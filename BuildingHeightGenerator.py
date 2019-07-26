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


args = parser.parse_args()

dtmFile = args.dtmFile
domFile = args.domFile
gridSquareSize = args.gridSquareSize
overpassRetryCount = args.overpassRetryCount
outputFilename = args.outputFilename
overpassFailRetryDelay = args.overpassFailRetryDelay

earthRadius = 6378.137
meterInDegree = ( 1 / ( ( 2 * math.pi / 360) * earthRadius ) ) / 1000

def downloadBuildings( left, bottom, right, top ):

    squareHeight = geopy.distance.geodesic( ( top, left ), ( bottom, left ) ).km
    squareWidth = geopy.distance.geodesic( ( top, left ), ( top, right ) ).km

    print(" \n\n\n Downloading buildings within the following square: ")
    print("  Bounds in epsg:4326: (left=%f, bottom=%f, right=%f, top=%f)" % ( ( left, bottom, right, top ) ) )
    print("  Bounds height: %f km, width: %f km" % ( squareHeight, squareWidth ) )

    api = overpy.Overpass()

    query = '''[out:json][timeout:500][maxsize:1000000000][bbox:%f,%f,%f,%f];
    ( 
      way["building"]; 
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

    bottom, right = sourceProj( lidar_dtm.bounds.left, lidar_dtm.bounds.top, inverse=True )
    top, left =  sourceProj( lidar_dtm.bounds.right, lidar_dtm.bounds.bottom, inverse=True )

    boundsHeight = geopy.distance.geodesic( ( left, top ), ( left, bottom ) ).km
    boundsWidth = geopy.distance.geodesic( ( left,top ), ( right,top ) ).km

    print("\n\n")
    print("Bounds in epsg:4326: (left=%f, bottom=%f, right=%f, top=%f)" % ( ( left, bottom, right, top ) ) )
    print("Bounds height: %f km, width: %f km" % ( boundsHeight, boundsWidth ) )


    dtmData = lidar_dtm.read(1)
    domData = lidar_dom.read(1)

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
    for x in range(longitudeGridIntervals-1):
        for y in range(latitudeGridIntervals-1):


            print("Starting x,y=" + str( ( x, y ) ) )


            gridCellLeft = addMeterLongitude( left, bottom, x * gridSquareSize * 1000 )
            gridCellTop = addMeterLatitude( left, bottom, ( y + 1 ) * gridSquareSize * 1000 )
            gridCellRight = addMeterLongitude( left, bottom, ( x + 1 ) * gridSquareSize * 1000 )
            gridCellBottom = addMeterLatitude( left, bottom, y * gridSquareSize * 1000 )
            
            buildingsResult = downloadBuildings( gridCellLeft, gridCellBottom, gridCellRight, gridCellTop )

            print("\nNumber of buildings downloaded: " + str(len( buildingsResult.ways )))
            for building in buildingsResult.ways:

                leftMostCol = math.inf
                rightMostCol = 0
                topMostRow = 0
                bottomMostRow = math.inf

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

                    nodesRowCol.append( ( row, col ) )

                    leftMostCol = min( leftMostCol, col )
                    rightMostCol = max( rightMostCol, col )
                    topMostRow = max( topMostRow, row )
                    bottomMostRow = min( bottomMostRow, row )

                averageLat = latitudeSum / len(building.nodes)
                averageLon = longitudeSum / len(building.nodes)

                #print(str(nodesRowCol))

                dtmSection = dtmData[ bottomMostRow:topMostRow, leftMostCol:rightMostCol ]
                domSection = domData[ bottomMostRow:topMostRow, leftMostCol:rightMostCol ]

                if(dtmSection.shape[0] == 0 or dtmSection.shape[1] == 0):
                    continue

                surfaceHeight = domSection - dtmSection

                buildingOutline = Image.new( "1", (  rightMostCol - leftMostCol, topMostRow - bottomMostRow  ), 0 )
                draw = ImageDraw.Draw(buildingOutline)

                pxs = list()
                for ( row, col ) in nodesRowCol:
                    pxs.append( ( col - leftMostCol, row - bottomMostRow ) )

                draw.polygon( pxs, 1, 1 )

                buildingOutlineArray = np.array( buildingOutline )


                #imgplot = plt.imshow(domSection)
                #plt.show()
                #imgplot = plt.imshow(dtmSection)
                #plt.show()

                #f, axarr = plt.subplots(1,4)
                #axarr[0].imshow(domSection)
                #axarr[1].imshow(dtmSection)
                #axarr[2].imshow(surfaceHeight)
                #axarr[3].imshow(buildingOutlineArray)

                #print(str(domSection - dtmSection))

                #plt.show()

                #input("Press Enter to continue...")

                #print(str(building.id))

                if(surfaceHeight[buildingOutlineArray].shape == (0,)):
                    continue

                buildingHeight = np.max( surfaceHeight[buildingOutlineArray] )
                wayUrl = "https://www.openstreetmap.org/way/" + str( building.id ) + ": " + str(buildingHeight) 
                print(wayUrl)


                if(buildingHeight > 0):

                    totalNumberOfBuildings = totalNumberOfBuildings + 1

                    outputFile.write('\n<node id="%d" lat="%f" lon="%f">' % ( -totalNumberOfBuildings, averageLat, averageLon ) )
                    outputFile.write('\n  <tag k="height" v="%f" />' % buildingHeight)
                    outputFile.write('\n</node>')

                    with open('buildings.txt', 'a') as f:
                        wayUrl = "https://www.openstreetmap.org/way/" + str( building.id ) + ": " + str(buildingHeight) 
                        f.write("\n" + wayUrl)


    outputFile.write('</osm>')


    print("\n\n\nTotal number of buildings: %d" % totalNumberOfBuildings )

