import rasterio
import os
import argparse
import overpy
from pyproj import Proj
from pyproj import transform
from functools import partial
import geopy.distance
import math
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageDraw
import time



parser = argparse.ArgumentParser()

parser.add_argument('dtmFile', type=str)
parser.add_argument('domFile', type=str)
parser.add_argument('outputFilename', type=str)
parser.add_argument('--gridSquareSize', type=float, default=2)
parser.add_argument('--overpassRetryCount', type=int, default=100)
parser.add_argument('--overpassFailRetryDelay', type=int, default=15)
parser.add_argument('--overpassApiUrl', type=str, default=None)


args = parser.parse_args()

dtmFile = args.dtmFile
domFile = args.domFile
gridSquareSize = args.gridSquareSize
overpassRetryCount = args.overpassRetryCount
outputFilename = args.outputFilename
overpassFailRetryDelay = args.overpassFailRetryDelay
overpassApiUrl = args.overpassApiUrl

earthRadius = 6378.137
meterInDegree = ( 1 / ( ( 2 * math.pi / 360) * earthRadius ) ) / 1000



def downloadBuildings( left, bottom, right, top ):

    squareHeight = geopy.distance.geodesic( ( top, left ), ( bottom, left ) ).km
    squareWidth = geopy.distance.geodesic( ( top, left ), ( top, right ) ).km

    print(" \n\n\n Downloading buildings within the following square: ")
    print("  Bounds in epsg:4326: (left=%f, bottom=%f, right=%f, top=%f)" % ( ( left, bottom, right, top ) ) )
    print("  Bounds height: %f km, width: %f km" % ( squareHeight, squareWidth ) )

    print("http://bboxfinder.com  bbox: ( %f, %f, %f, %f )" % ( left, bottom, right, top ))

    api = overpy.Overpass(url=overpassApiUrl)

    query = '''[out:json][timeout:500][maxsize:1000000000][bbox:%f,%f,%f,%f];
    ( 
      way["building"]; 
      way["building:part"]; 
    );
    (._;>;);
    out meta;''' % ( bottom, left, top, right )

    
    #print(query)

    for i in range(overpassRetryCount):
        try:
            
            result = api.query(query)

            return result
        except:
            print("Failed, retrying in %d seconds" % overpassFailRetryDelay)
            time.sleep(overpassFailRetryDelay)
    
    print("Unable to query overpass with %d retries" % overpassRetryCount)
    exit()
        


with rasterio.open(dtmFile) as lidar_dtm,  rasterio.open(domFile) as lidar_dom, open(outputFilename, "w") as outputFile:

    if lidar_dtm.bounds != lidar_dom.bounds:
        print("DTM and DOM does not cover the same extent")
        exit()

    sourceProj = Proj( init='EPSG:25833' )
    wgsProj = Proj( init='EPSG:4326' )
    
    sourceTransformer = partial(transform, wgsProj, sourceProj)

    print("Bounds in source coordinate system" + str(lidar_dtm.bounds))

    ( leftTopCornerLong, leftTopCornerLat ) = sourceProj( lidar_dtm.bounds.left, lidar_dtm.bounds.top, inverse=True )
    ( rightTopCornerLong, rightTopCornerLat ) = sourceProj( lidar_dtm.bounds.right, lidar_dtm.bounds.top, inverse=True )
    ( leftBottomCornerLong, leftBottomCornerLat ) = sourceProj( lidar_dtm.bounds.left, lidar_dtm.bounds.bottom, inverse=True )
    ( rightBottomCornerLong, rightBottomCornerLat ) = sourceProj( lidar_dtm.bounds.right, lidar_dtm.bounds.bottom, inverse=True )

    left = min( leftTopCornerLong, leftBottomCornerLong )
    right = max( rightTopCornerLong, rightBottomCornerLong )
    bottom = min( leftBottomCornerLat,rightBottomCornerLat )
    top = max( leftTopCornerLat, rightTopCornerLat )


    boundsHeight = geopy.distance.geodesic( ( top, left ), ( bottom, left ) ).km
    boundsWidth = geopy.distance.geodesic( ( top, left ), ( top, right ) ).km

    print("\n\n")
    print("Bounds in epsg:4326: (left=%f, bottom=%f, right=%f, top=%f)" % ( ( left, bottom, right, top ) ) )
    print("Bounds height: %f km, width: %f km" % ( boundsHeight, boundsWidth ) )
    print("http://bboxfinder.com  bbox: ( %f, %f, %f, %f )" % ( left, bottom, right, top ))

    print("Starting to load DTM")
    dtmData = lidar_dtm.read(1)
    print("DTM loaded")

    print("Starting to load DSM")
    domData = lidar_dom.read(1)
    print("DSM Loaded")

    #print(str(lidar_dom.meta))


    longitudeGridIntervals = math.ceil(  boundsWidth / gridSquareSize )
    latitudeGridIntervals = math.ceil( boundsHeight / gridSquareSize )

    print( "Latitude intervals: %d, longitude intervals: %d " % ( latitudeGridIntervals, longitudeGridIntervals ) )
    print( "Total number of squares: %d" % ( latitudeGridIntervals * longitudeGridIntervals )  )

    
    outputFile.write('<?xml version="1.0" encoding="UTF-8"?>')
    outputFile.write('\n<osm version="0.6" generator="BuildingHeightGenerator">')

    def addMeterLongitude( latitude, longitude, meters):

        return longitude + (meters * meterInDegree) / math.cos( latitude * ( math.pi / 180 ) )

    def addMeterLatitude( latitude, longitude, meters ):
        
        return latitude + (meters * meterInDegree)

    totalNumberOfBuildings = 0
    for x in range(longitudeGridIntervals):
        for y in range(latitudeGridIntervals):


            print("Starting x,y=" + str( ( x, y ) ) )


            gridCellLeft = addMeterLongitude( bottom, left, x * gridSquareSize * 1000 )
            gridCellTop = addMeterLatitude( bottom, left, ( y + 1 ) * gridSquareSize * 1000 )
            gridCellRight = addMeterLongitude( bottom, left, ( x + 1 ) * gridSquareSize * 1000 )
            gridCellBottom = addMeterLatitude( bottom, left, y * gridSquareSize * 1000 )

            #print("http://bboxfinder.com  grid_bbox: ( %f, %f, %f, %f )" % ( gridCellLeft, gridCellBottom, gridCellRight, gridCellTop ))


            def checkNoDataInGridCell( data, noDataValue ):

                ( gridCellLeftX, gridCellTopY ) = sourceTransformer( gridCellLeft, gridCellTop )
                ( gridCellRightX, gridCellBottomY )  = sourceTransformer( gridCellRight, gridCellBottom )


                ( gridCellTopRow, gridCellLeftCol ) = lidar_dtm.index( gridCellLeftX, gridCellTopY )
                ( gridCellBottomRow, gridCellRightCol ) = lidar_dtm.index( gridCellRightX, gridCellBottomY )

                gridCellData = data[ max(0,gridCellTopRow):max(0,gridCellBottomRow), max(0, gridCellLeftCol):max(0,gridCellRightCol) ]

                return np.all( gridCellData == noDataValue )

            
            if( checkNoDataInGridCell(dtmData, lidar_dtm.nodata) ):
                print("grid contains only no-data in DTM file: skipping")
                continue

            if( checkNoDataInGridCell(domData, lidar_dom.nodata) ):
                print("grid contains only no-data in DOM file: skipping")
                continue
            
            
            buildingsResult = downloadBuildings( gridCellLeft, gridCellBottom, gridCellRight, gridCellTop )

            print("\nNumber of buildings downloaded: " + str(len( buildingsResult.ways )))
            for building in buildingsResult.ways:

                leftMostCol = math.inf
                rightMostCol = 0
                topMostRow = math.inf
                bottomMostRow = 0

                nodesRowCol = []
                latitudeSum = 0
                longitudeSum = 0
                for node in building.nodes:

                    latitudeSum = latitudeSum + node.lat 
                    longitudeSum = longitudeSum + node.lon
                    
                    #print(str(node))
                    coordinate = sourceTransformer( node.lon, node.lat )
                    xCoordinate,yCoordinate = coordinate

                    #print(str(coordinate))

                    row, col = lidar_dtm.index( xCoordinate, yCoordinate )

                    #print(str((row, col)))
                    row = max(0, row)
                    col = max(0, col)

                    nodesRowCol.append( ( row, col ) )

                    leftMostCol = min( leftMostCol, col )
                    rightMostCol = max( rightMostCol, col )
                    topMostRow = min( topMostRow, row )
                    bottomMostRow = max( bottomMostRow, row )

                averageLat = latitudeSum / len(building.nodes)
                averageLon = longitudeSum / len(building.nodes)

                #print(str(nodesRowCol))

                dtmSection = dtmData[ topMostRow:bottomMostRow, leftMostCol:rightMostCol ]
                domSection = domData[ topMostRow:bottomMostRow, leftMostCol:rightMostCol ]


                ( sectionHeight, sectionWidth ) = dtmSection.shape

                if(sectionHeight == 0 or sectionWidth == 0):
                    continue

                surfaceHeight = domSection - dtmSection

                buildingOutline = Image.new( "1", (  sectionWidth, sectionHeight  ), 0 )
                draw = ImageDraw.Draw(buildingOutline)

                pxs = list()
                for ( row, col ) in nodesRowCol:
                    pxs.append( ( col - leftMostCol, row - topMostRow ) )

                draw.polygon( pxs, 1, 1 )

                buildingOutlineArray = np.array( buildingOutline )


                #f, axarr = plt.subplots(1,4)
                #axarr[0].imshow(domSection)
                #axarr[1].imshow(dtmSection)
                #axarr[2].imshow(surfaceHeight)
                #axarr[3].imshow(buildingOutlineArray)

                #plt.show()

                #input("Press Enter to continue...")
                #exit()

                #print(str(building.id))

                buildingOutlineSurfaceHeight = surfaceHeight[buildingOutlineArray]

                if(buildingOutlineSurfaceHeight.size == 0):
                    continue

                buildingHeight = np.max( buildingOutlineSurfaceHeight )
                wayUrl = "https://www.openstreetmap.org/way/" + str( building.id ) + ": " + str(buildingHeight) 
                print(wayUrl)

                
                if(buildingHeight > 1):

                    buildingHeight = round(buildingHeight,1)

                    totalNumberOfBuildings = totalNumberOfBuildings + 1

                    outputFile.write('\n<node id="%d" lat="%f" lon="%f">' % ( -totalNumberOfBuildings, averageLat, averageLon ) )
                    outputFile.write('\n  <tag k="height" v="%.1f" />' % buildingHeight)
                    outputFile.write('\n</node>')

                    with open('buildings.txt', 'a') as f:
                        wayUrl = "https://www.openstreetmap.org/way/" + str( building.id ) + ": " + str(buildingHeight) 
                        f.write("\n" + wayUrl)

    outputFile.write('</osm>')


    print("\n\n\nTotal number of buildings: %d" % totalNumberOfBuildings )

