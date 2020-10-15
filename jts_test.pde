import org.locationtech.jts.io.*;
import org.locationtech.jts.geom.*;
import org.locationtech.jts.algorithm.*;
import org.locationtech.jts.geom.util.*;
import org.locationtech.jts.operation.distance.*;
import org.locationtech.jts.noding.*;
import org.locationtech.jts.algorithm.*;
import org.locationtech.jts.geom.prep.*;
import org.locationtech.jts.triangulate.*;
import java.util.*;
import extruder.*;


Geometry g1;
Polygon outline;

Coordinate firstPoint = new Coordinate(0, 0);
Coordinate lastPoint = new Coordinate(0, 0);

int f = 0; //frames counter

void setup() {
  size(800, 800, P3D);
  smooth();


  try {

    //Coordinate[] coords  = new Coordinate[] {new Coordinate(400, 0),  new Coordinate(200, 200),  new Coordinate(400, 400), new Coordinate(600, 200), new Coordinate(400, 0) };
    //Polygon polygon = new GeometryFactory().createPolygon(coords);

    g1 = new WKTReader().read("LINESTRING (-77 51, 10 130, 40 140, 60 130, 37 126)");

    outline = pointToCircle( new Coordinate(0, 0), 300 );
    //outline = (Polygon)new WKTReader().read("POLYGON ((-100 100, 100 100, 0 -100, -100 100))");
  } 
  catch (Exception ex) {
  }
}

Polygon[] splitOutline(LineString curve, LineString outlineLineString, GeometryLocation firstLocation, GeometryLocation lastLocation) {
  Coordinate[] outlineCoords = outlineLineString.getCoordinates();

  List<Coordinate> coords1 = new ArrayList();
  coords1.addAll((List<Coordinate>)Arrays.asList(curve.getCoordinates()));
  coords1.add(lastLocation.getCoordinate());

  int segmentIndex1 = lastLocation.getSegmentIndex();
  while (segmentIndex1!=firstLocation.getSegmentIndex()) {
    coords1.add(outlineCoords[segmentIndex1+1]);
    segmentIndex1 = (segmentIndex1 + 1) % (outlineCoords.length-1);
  }


  coords1.add(firstLocation.getCoordinate());
  coords1.add(curve.getCoordinates()[0]);  
  Polygon polygon1 = (new GeometryFactory()).createPolygon(coords1.toArray(new Coordinate[0]));


  List<Coordinate> coords2 = new ArrayList();
  List<Coordinate> curveCoords = (List<Coordinate>)Arrays.asList(curve.getCoordinates());
  Collections.reverse(curveCoords);
  coords2.addAll(curveCoords);
  coords2.add(firstLocation.getCoordinate());

  int segmentIndex2 = firstLocation.getSegmentIndex();
  while (true) {
    coords2.add(outlineCoords[segmentIndex2+1]);
    segmentIndex2 = (segmentIndex2 + 1) % (outlineCoords.length-1);
    if (segmentIndex2==lastLocation.getSegmentIndex()) {
      break;
    }
  }

  coords2.add(lastLocation.getCoordinate());
  coords2.add(curveCoords.get(0));  
  Polygon polygon2 = (new GeometryFactory()).createPolygon(coords2.toArray(new Coordinate[0]));

  // take buffer to simplify shapes -- can return multiple polygons
  //GeometryCollection geometries1 = new GeometryCollection(new Geometry[] { polygon1 }, new GeometryFactory());
  //polygon1 = (Polygon)geometries1.buffer(0);

  //GeometryCollection geometries2 = new GeometryCollection(new Geometry[] { polygon2 }, new GeometryFactory());
  //polygon2 = (Polygon)geometries2.buffer(0);

  return new Polygon[] {
    polygon1, polygon2
  };
}

Geometry bufferPolygon(Polygon p, float distance) {
  GeometryCollection geometries2 = new GeometryCollection(new Geometry[] { p }, new GeometryFactory());
  return geometries2.buffer(distance);
}

PShape geometryToPShape(Geometry g) {
  Coordinate[] coords = g.getCoordinates();
  PShape s = createShape();
  s.beginShape();
  for (int i=0; i<coords.length; i++) {
    s.vertex((float)coords[i].x, (float)coords[i].y);
  }
  s.endShape();
  return s;
}

Polygon pointToCircle(Coordinate coord, float diameter) {
  GeometricShapeFactory factory = new GeometricShapeFactory();
  factory.setCentre(coord);
  factory.setWidth(diameter);
  factory.setHeight(diameter);
  return factory.createCircle();
}

Geometry coordToGeometry(Coordinate coord) {
  return (new GeometryFactory()).createPoint(coord);
}

LineString unringGeometry(Geometry geom) {
  Coordinate[] srcCoords = geom.getCoordinates();
  return new GeometryFactory().createLineString(srcCoords);
}

void draw() {
  background(255);

  noFill();
  stroke(255);

  if (mousePressed) {
    if (mouseButton == LEFT) {
      firstPoint = new Coordinate((double)mouseX - width/2, (double)mouseY - height/2);
    } else {
      lastPoint = new Coordinate((double)mouseX - width/2, (double)mouseY - height/2);
    }
  }

  Coordinate[] coordinates = new Coordinate[g1.getCoordinates().length+2];
  System.arraycopy(g1.getCoordinates(), 0, coordinates, 1, g1.getCoordinates().length);

  coordinates[0] = firstPoint;
  coordinates[coordinates.length-1] = lastPoint;

  LineString curve = new GeometryFactory().createLineString(coordinates);

  pushMatrix();
  translate(width/2, height/2);







  MinimumBoundingCircle mb = new MinimumBoundingCircle(curve);
  Geometry boundingCircle = mb.getCircle();
  Coordinate center = mb.getCentre();
  boundingCircle = (new AffineTransformation().translate(-center.x, -center.y).scale(1.2, 1.2).translate(center.x, center.y)).transform(boundingCircle);


  shape(geometryToPShape(curve));
  stroke(255, 255, 255);
  shape(geometryToPShape(outline));


  pushStyle();
  stroke(127);
  shape(geometryToPShape(boundingCircle));
  popStyle();


  boolean firstPointInside = outline.contains(coordToGeometry(firstPoint));
  boolean lastPointInside = outline.contains(coordToGeometry(lastPoint));

  pushStyle();
  fill(firstPointInside ? 255 : 0, 0, 0 );
  shape(geometryToPShape(
    pointToCircle(firstPoint, 10)
    ));

  fill(lastPointInside ? 255 : 0, 0, 0 );
  shape(geometryToPShape(
    pointToCircle(lastPoint, 10)
    ));
  popStyle();

  // first connecting points on the outline
  Coordinate firstConnect = new Coordinate(0, 0);
  Coordinate lastConnect = new Coordinate(0, 0);


  List<SegmentString> curveStrings = (List<SegmentString>)SegmentStringUtil.extractSegmentStrings(curve);

  LineString outlineLineString = unringGeometry(outline);
  //if (firstPointInside) {
  // find nearest point on the outline
  DistanceOp distanceOp = new DistanceOp(coordToGeometry(firstPoint), outlineLineString);
  GeometryLocation[] firstLocations = distanceOp.nearestLocations();
  GeometryLocation firstLocation = firstLocations[1]; 
  firstConnect = firstLocation.getCoordinate();
  //} else {

  //  PreparedLineString preparedOutline = new PreparedLineString(unringGeometry(outline));
  //  SegmentIntersectionDetector detector = new SegmentIntersectionDetector();

  //  preparedOutline.getIntersectionFinder().intersects(curveStrings, detector);
  //  firstConnect = detector.getIntersection();   

  //}


  //if (lastPointInside) {
  // find nearest point on the outline
  DistanceOp lastDistanceOp = new DistanceOp(coordToGeometry(lastPoint), outlineLineString);

  GeometryLocation[] lastLocations = lastDistanceOp.nearestLocations();
  GeometryLocation lastLocation = lastLocations[1]; 
  lastConnect = lastLocation.getCoordinate();
  println("first: " + firstLocation.getSegmentIndex() + " last: " + lastLocation.getSegmentIndex());
  //} else {

  //  PreparedLineString preparedOutline = new PreparedLineString(unringGeometry(outline));
  //  SegmentIntersectionDetector detector = new SegmentIntersectionDetector();

  //  preparedOutline.getIntersectionFinder().intersects(curveStrings, detector);
  //  lastConnect = detector.getIntersection();    
  //}

  //Geometry intersectionGeom = outline.difference(curve);

  //stroke(0,100,100);
  //fill(255,0,0);
  //shape(geometryToPShape(intersectionGeom));



  pushStyle();
  fill(0, 127, 0);
  shape(geometryToPShape(
    pointToCircle(firstConnect, 10)
    ));

  fill(0, 0, 127);
  shape(geometryToPShape(
    pointToCircle(lastConnect, 10)
    ));
  popStyle();

  IntersectionMatrix matrix = curve.relate(outline);
  text(matrix.toString(), 10, -100);

  Polygon[] splitPolygons = splitOutline(curve, outlineLineString, firstLocation, lastLocation);

  Geometry splitGeometry1 = bufferPolygon(splitPolygons[0], 0).intersection(outline);

  Geometry splitGeometry2 = bufferPolygon(splitPolygons[1], 0).intersection(outline);


  stroke(0, 100, 100);
  noFill();
  shape(geometryToPShape(splitPolygons[0]));
  fill(255, 0, 0, 127);
  shape(geometryToPShape(splitGeometry1));
  noFill();
  //shape(geometryToPShape(splitPolygons[1]));
  fill(0, 255, 0, 127);
  //shape(geometryToPShape(splitGeometry2));

  //GeometryFactory geomFact = new GeometryFactory();
  //DelaunayTriangulationBuilder builder = new DelaunayTriangulationBuilder();
  //builder.setSites(splitGeometry1);

  //Geometry result = null;
  ////if (true) {
  ////result = builder.getTriangles(geomFact);
  ////}
  ////else {
  //  result = builder.getEdges(geomFact);
  ////}
  PShape[] extruded;

  fill(100, 255, 100, 127);

  extruder e;
  e = new extruder(this);
  extruded = e.extrude(geometryToPShape(splitPolygons[0]), 140, "box");

  //shape(geometryToPShape(splitPolygons[0]));

  //shape(geometryToPShape(result));

  //System.out.println(result);
  //background(0);

  // Set origin of scene to center of image
  // Rotate 3 degrees per frame on the y-axis
  rotateY(radians(f*3));

  PShape _shape= createShape(GROUP);
  for (PShape p: extruded){
    _shape.addChild(p);
  }
  shape(_shape);
  // Increment frame counter
  f++;

  popMatrix();
}
