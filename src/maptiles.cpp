/**
 * @file maptiles.cpp
 * Code for the maptiles function.
 */

#include <iostream>
#include <map>

#include "maptiles.h"

using namespace std;


Point<3> convertToXYZ(LUVAPixel pixel) {
    return Point<3>( pixel.l, pixel.u, pixel.v );
}

MosaicCanvas* mapTiles(SourceImage const& theSource,
                       vector<TileImage>& theTiles)
{
    /**
     * @todo Implement this function!
     */
    vector<Point<3>> points;
    map<Point<3>, TileImage*> tile_map;
    for (unsigned i = 0; i < theTiles.size(); i++) {
        points.push_back(convertToXYZ(theTiles[i].getAverageColor()));
        tile_map[convertToXYZ(theTiles[i].getAverageColor())] = &theTiles[i];
    }

    KDTree<3> tree(points);

    MosaicCanvas * canvas = new MosaicCanvas(theSource.getRows(), theSource.getColumns());
    for (int r = 0; r < theSource.getRows(); r++) {
        for (int c = 0; c < theSource.getColumns(); c++) {
            LUVAPixel original = theSource.getRegionColor(r, c);
            Point<3> neighbor = tree.findNearestNeighbor(convertToXYZ(original));
            TileImage * tile = tile_map[neighbor];
            canvas->setTile(r, c, tile);
        }
    }
    return canvas;
}

