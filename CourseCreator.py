import math, pyproj, numpy
from xml.etree.ElementTree import Element, SubElement, Comment, tostring
import xml.etree.ElementTree as ET
import random
import requests
import zipfile
import tempfile
from shapely import geometry
from wand.image import Image,Color

# Server for route distance calculation
OSRM_SERVER = "192.168.1.85"
# Server for map creation
PDF_SERVER = "http://192.168.1.85:8080"

KML_NAMESPACE = "http://www.opengis.net/kml/2.2"

class Point:
    '''A location on the earth.'''
    def __init__(self, lat, lon):
        self.lat = float(lat)
        self.lon = float(lon)

    def __eq__(self, other):
        return self.lat == other.lat and self.lon == other.lon

    def __hash__(self):
        return hash((self.lat,self.lon))

    def getLat(self):
        return self.lat

    def getLon(self):
        return self.lon

    def getEPSG3857(self):
        '''Get coordinates in the format required to request a map.'''
        inProj = pyproj.Proj(init='epsg:4326')
        outProj = pyproj.Proj(init='epsg:3857')
        x, y = pyproj.transform(inProj, outProj, self.lon, self.lat)
        return (x, y)

    def getOffsetFromPoint(self, point2):
        '''Calculate the (x,y) offset in metres from another point.'''
        fudgeFactor = math.cos(point2.getLat() * math.pi / 180)
        x1, y1 = self.getEPSG3857()
        x2, y2 = point2.getEPSG3857()
        return (fudgeFactor * (x1 - x2), fudgeFactor * (y1 - y2))

    @classmethod
    def fromEPSG3857(cls, x, y):
        '''Create from EPSG3857 coordinates.'''
        outProj = pyproj.Proj(init='epsg:4326')
        inProj = pyproj.Proj(init='epsg:3857')
        lon, lat = pyproj.transform(inProj, outProj, x, y)
        return cls(lat, lon)

    @classmethod
    def fromBNG(cls, easting, northing):
        '''Create from British National Grid.'''
        outProj = pyproj.Proj(init='epsg:4326')
        inProj = pyproj.Proj(init='epsg:27700')
        lon, lat = pyproj.transform(inProj, outProj, easting, northing)
        return cls(lat, lon)

class Control(Point):
    '''A Point with a control number.'''
    def __init__(self, point, number):
        super().__init__(point.lat, point.lon)
        self.number = number

    def getNumber(self):
        return self.number

    def getScore(self):
        # Use NGOC scoring system
        return 10 * math.ceil(self.number / 10)

class Course:
    '''The underlying Controls and start/finish.'''
    def __init__(self, startFinish):
        self.startFinish = startFinish
        self.controlSet = set()

    def addControls(self, controlsToAdd):
        for control in controlsToAdd:
            self.controlSet.add(control)

    def getControls(self):
        return self.controlSet

class CourseRepresentation(Course):
    '''Everything needed to create the right files from a Course.'''
    def __init__(self, course):
        self.startFinish = course.startFinish
        self.controlSet = course.controlSet
        print([control for control in course.controlSet])
        xCoords = [control.getEPSG3857()[0] for control in course.controlSet]
        yCoords = [control.getEPSG3857()[1] for control in course.controlSet]
        midX = (max(xCoords) + min(xCoords)) / 2
        midY = (max(yCoords) + min(yCoords)) / 2
        spreadX = max(xCoords) - min(xCoords)
        spreadY = max(yCoords) - min(yCoords)
        self.midpoint = Point.fromEPSG3857(midX, midY)
        print ("Spread",spreadX,spreadY)
        if (spreadX > spreadY):
            self.xDim = 297
            self.yDim = 210
            self.aspect = "landscape"
        else:
            self.xDim = 210
            self.yDim = 297
            self.aspect = "portrait"
        fudgeFactor = math.cos(self.midpoint.getLat() * math.pi / 180)
        self.topLeft = Point.fromEPSG3857(midX-(10.0/fudgeFactor)*(self.xDim/2),midY + (10.0/fudgeFactor)*(self.yDim/2))
        self.bottomRight = Point.fromEPSG3857(midX+(10.0/fudgeFactor)*(self.xDim/2),midY - (10.0/fudgeFactor)*(self.yDim/2))

    def makePpenFile(self):
        '''Create PurplePen file, for nicer printing. Apparently this might be produced by the MapRun server anyway.'''
        courseScribeEvent = Element("course-scribe-event")
        event = SubElement(courseScribeEvent, "event", {"id": "1"})
        #TODO Get the pdf
        map = SubElement(event, "map", {"kind": "PDF", "scale": "10000", "absolute-path": "PATH"})
        allControls = SubElement(event, "all-controls", {"description-kind": "symbols", "print-scale": "10000"})
        if self.aspect == "landscape":
            printArea = SubElement(event, "print-area",
                                   {"automatic": "true", "restrict-to-page-size": "true", "left": "0", "top": "210",
                                    "right": "297", "bottom": "0", "page-width": "827", "page-height": "1169",
                                    "page-margins": "0", "page-landscape": "true"})
        else:
            printArea = SubElement(event, "print-area",
                                   {"automatic": "true", "restrict-to-page-size": "true", "left": "0", "top": "297",
                                    "right": "210", "bottom": "0", "page-width": "827", "page-height": "1169",
                                    "page-margins": "0", "page-landscape": "false"})
        numbering = SubElement(event, "numbering", {"start": "1", "disallow-invertible": "false"})
        courseAppearance = SubElement(event, "course-appearance",
                                      {"center-dot-diameter": "1", "auto-leg-gap-size": "3.5", "blend-purple": "true"})
        startFinishX = self.xDim/2 + self.startFinish.getOffsetFromPoint(self.midpoint)[0]/10
        startFinishY = self.yDim/2 + self.startFinish.getOffsetFromPoint(self.midpoint)[1]/10
        print (self.startFinish.getOffsetFromPoint(self.midpoint)[1])
        startControl = SubElement(courseScribeEvent, "control", {"id": "101", "kind": "start"})
        startLocation = SubElement(startControl, "location", {"x":str(startFinishX),"y":str(startFinishY)})
        finishControl = SubElement(courseScribeEvent, "control", {"id": "102", "kind": "finish"})
        finishDescription = SubElement(finishControl,"description",{"box":"all", "iof-2004-ref":"14.3"})
        finishLocation = SubElement(finishControl, "location", {"x":str(startFinishX),"y":str(startFinishY)})
        i=1
        normalControls=dict()
        for controlObject in self.controlSet:
            normalControls[i] = SubElement(courseScribeEvent, "control", {"id": str(i), "kind": "normal"})
            code = SubElement(normalControls[i],"code")
            code.text = str(controlObject.getNumber())
            normalX = self.xDim/2 + controlObject.getOffsetFromPoint(self.midpoint)[0]/10
            normalY = self.yDim/2 + controlObject.getOffsetFromPoint(self.midpoint)[1]/10
            normalLocation = SubElement(normalControls[i], "location", {"x": str(normalX), "y": str(normalY)})
            i+=1
        courseObject=SubElement(courseScribeEvent,"course", {"id":"1","kind":"score","order":"1"})
        courseName=SubElement(courseObject,"name")
        courseName.text="Name"
        labels = SubElement(courseObject,"labels",{"label-kind":"code"})
        first = SubElement(courseObject,"first",{"course-control":"101"})
        if self.aspect == "landscape":
            printArea = SubElement(courseObject, "print-area",
                                   {"automatic": "true", "restrict-to-page-size": "true", "left": "0", "top": "210",
                                    "right": "297", "bottom": "0", "page-width": "827", "page-height": "1169",
                                    "page-margins": "0", "page-landscape": "true"})
        else:
            printArea = SubElement(courseObject, "print-area",
                                   {"automatic": "true", "restrict-to-page-size": "true", "left": "0", "top": "297",
                                    "right": "210", "bottom": "0", "page-width": "827", "page-height": "1169",
                                    "page-margins": "0", "page-landscape": "false"})
        startObject = SubElement(courseScribeEvent,"course-control",{"id":"101","control":"101"})
        startNext = SubElement(startObject,"next",{"course-control":"1"})
        finishObject = SubElement(courseScribeEvent,"course-control",{"id":"102","control":"102"})
        controlObjects=dict()
        i = 1
        for controlObject in self.controlSet:
            controlObjects[i]=SubElement(courseScribeEvent,"course-control", {"id":str(i),"control":str(i)})
            if i==len(self.controlSet):
                nextControl = SubElement(controlObjects[i],"next",{"course-control":"102"})
            else:
                nextControl = SubElement(controlObjects[i],"next",{"course-control":str(i+1)})
            i+=1
        print(tostring(courseScribeEvent))
        print(self.midpoint.getLat(), self.midpoint.getLon())

    def makeKmlFile(self,filename):
        '''Produce the course for MapRun: KML format.'''
        kml = Element("kml", {"xmlns":"http://www.opengis.net/kml/2.2"})
        document = SubElement(kml,"Document")
        folder = SubElement(document,"Folder")
        startPlacemark=SubElement(folder,"Placemark")
        name=SubElement(startPlacemark,"name")
        name.text="S1"
        point = SubElement(startPlacemark,"Point")
        coordinates=SubElement(point,"coordinates")
        coordinates.text="%f,%f,0" % (self.startFinish.getLon(),self.startFinish.getLat())
        placemarks = dict()
        for controlObject in self.controlSet:
            placemarks[controlObject.getNumber()]=SubElement(folder,"Placemark")
            name=SubElement(placemarks[controlObject.getNumber()],"name")
            name.text=str(controlObject.getNumber())
            point = SubElement(placemarks[controlObject.getNumber()],"Point")
            coordinates=SubElement(point,"coordinates")
            coordinates.text="%f,%f,0" % (controlObject.getLon(),controlObject.getLat())
            controlObject.getNumber()
        finishPlacemark=SubElement(folder,"Placemark")
        name=SubElement(finishPlacemark,"name")
        name.text="F1"
        point = SubElement(finishPlacemark,"Point")
        coordinates=SubElement(point,"coordinates")
        coordinates.text="%f,%f,0" % (self.startFinish.getLon(),self.startFinish.getLat())
        f=open(filename, "wb")
        f.write(tostring(kml))
        f.close()

    def makeKmzFile(self,filename,title):
        '''Produce the MapRun map, as a georeferenced JPG image in a KMZ file.'''
        if self.aspect == "landscape":
            ss = "style=streeto|paper=0.297,0.210|scale=10000|centre="
        else:
            ss = "style=streeto|paper=0.210,0.297|scale=10000|centre="
        mp = self.getMidpoint().getEPSG3857()
        ss += str(int(mp[1])) + "," + str(int(mp[0])) + "|title="+title+"|club=|mapid=|start="
        ss += str(int(self.startFinish.getEPSG3857()[1])) + "," + str(int(self.startFinish.getEPSG3857()[0]))
        ss += "|crosses=|controls="
        #TODO Remove controls from request for actual MapRun map
        ss += ",".join([str(control.getNumber()) + ",45," + str(int(control.getEPSG3857()[1])) + "," + str(
            int(control.getEPSG3857()[0])) for control in self.getControls()])
        r = requests.get(PDF_SERVER+"/pdf/?"+ss)
        pdfFile = tempfile.NamedTemporaryFile(delete=False)
        for chunk in r.iter_content(chunk_size=128):
            pdfFile.write(chunk)
        pdfFile.close()
        with Image(filename=pdfFile.name, resolution=180) as img_pdf:
            img = Image(image=img_pdf.sequence[0])
            img.background_color = Color('white')
            img.alpha_channel = 'remove'
            img.format = "jpeg"
            img.save(filename=pdfFile.name+".jpg")
        kml = Element("kml", {"xmlns": "http://www.opengis.net/kml/2.2"})
        document = SubElement(kml, "Document")
        folder = SubElement(document, "Folder")
        groundOverlay = SubElement(folder, "GroundOverlay")
        name=SubElement(groundOverlay,"name")
        name.text="tile.jpg"
        drawOrder=SubElement(groundOverlay,"drawOrder")
        drawOrder.text="75"
        icon=SubElement(groundOverlay,"Icon")
        href=SubElement(icon,"href")
        href.text="files/tile.jpg"
        latLonBox=SubElement(groundOverlay,"LatLonBox")
        north=SubElement(latLonBox,"north")
        north.text=str(self.topLeft.getLat())
        south=SubElement(latLonBox,"south")
        south.text=str(self.bottomRight.getLat())
        east=SubElement(latLonBox,"east")
        east.text=str(self.bottomRight.getLon())
        west=SubElement(latLonBox,"west")
        west.text=str(self.topLeft.getLon())
        rotation=SubElement(latLonBox,"rotation")
        rotation.text="0.0"
        docKmlFile = tempfile.NamedTemporaryFile(delete=False)
        docKmlFile.write(tostring(kml))
        docKmlFile.close()
        with zipfile.ZipFile(filename,"w") as kmzFile:
            kmzFile.write(docKmlFile.name,"doc.kml")
            kmzFile.write(pdfFile.name+".jpg","files/tile.jpg")

    def getMidpoint(self):
        return self.midpoint

class Routes:
    '''Deals with the route responses from the OSRM server.'''
    def __init__(self,routesResponse):
        self.distances = []
        self.twists = []
        for route in routesResponse:
            self.distances.append(route['distance'])
            twistiness = 0
            for step in route['legs'][0]['steps']:
                maneuver = step['maneuver']
                if maneuver['type'] == 'turn':
                    twist = (maneuver['bearing_before']-maneuver['bearing_after'])%360
                    if twist > 180:
                        twist = 360 - twist
                    twistiness += twist
            self.twists.append(twistiness)

    def getDistances(self):
        '''Distances of the calculated routes (increasing).'''
        return self.distances

    def getTwists(self):
        '''A measure of how many turns are needed to execute the route.'''
        return self.twists

class Osrm:
    '''For Open Source Routing Machine server connections.'''
    def __init__(self, url):
        self.server = url

    def getNearest(self, point):
        '''Nearest waypoint (on the road network) to a point.'''
        r = requests.get(
            "http://"+self.server+":5000/nearest/v1/driving/" + str(point.getLon()) + "," + str(point.getLat()))
        return r.json()['waypoints'][0]

    def getRoutes(self, point1, point2):
        '''Routes between two points.  Returns an alternative if one found.'''
        r = requests.get(
            "http://"+self.server+":5000/route/v1/driving/" + str(point1.getLon()) + "," + str(point1.getLat()) + ";" +
            str(point2.getLon()) + "," + str(point2.getLat()), params={"steps": "true", "alternatives": "true"})
        return Routes(r.json()['routes'])

    def getRoundtripLength(self, course):
        '''Shortest distance to visit all points and return to the start.'''
        requestString = "http://"+self.server+":5000/trip/v1/driving/" + str(course.startFinish.getLon()) + "," + str(course.startFinish.getLat()) + ";"
        for control in course.getControls():
            requestString += str(control.getLon()) + "," + str(control.getLat()) + ";"
        requestString += str(course.startFinish.getLon()) + "," + str(course.startFinish.getLat())
        r = requests.get(requestString,params={"source":"first","destination":"last"})
        return r.json()['trips'][0]['distance']

class Candidates:
    '''The set of all possible points, from which the course is created.'''
    def __init__(self, startFinish):
        self.points = set()
        self.startFinish = startFinish

    def __init__(self, points, startFinish):
        self.points = set(points)
        self.startFinish = startFinish

    def addPoints(self,pointsToAdd):
        self.points.add(pointsToAdd)

    def courseByMinDistance(self, nControls, minDistance):
        '''Randomly select points > minDistance apart.'''
        osrm = Osrm(OSRM_SERVER)
        controls = []
        pointSet = set([self.startFinish])
        i=1
        while i<=nControls:
            newPoint = random.sample(self.points,1)[0]
            if osrm.getNearest(newPoint)['distance']>5:
                continue
            suitable = True
            for c in pointSet:
                if osrm.getRoutes(newPoint,c).getDistances()[0] < minDistance:
                    suitable = False
                    break
            if suitable:
                pointSet.add(newPoint)
                controls.append(Control(newPoint,i))
                i+=1
        course = Course(self.startFinish)
        course.addControls(controls)
        return course

    def courseByTwistiness(self, nControls):
        '''Pairs of points with route choice. Not a proper course!'''
        osrm = Osrm(OSRM_SERVER)
        controls = []
        i=1
        while i<=nControls:
            newPoints = random.sample(self.points,2)
            routes = osrm.getRoutes(newPoints[0],newPoints[1])
            twists = routes.getTwists()
            if len(twists)>1 and twists[0]>3*twists[1] and routes.getDistances()[0]<600:
                print(i)
                controls.append(Control(newPoints[0],i))
                controls.append(Control(newPoints[1],i+1))
                i+=2
        course = Course(self.startFinish)
        course.addControls(controls)
        return course

    def courseByAddition(self,nControls,minDistance,maxDistance):
        '''Add controls that are fairly equidistant from lots of others.'''
        controlOrder = numpy.random.permutation(range(1,31))
        osrm = Osrm(OSRM_SERVER)
        controls = []
        i=0
        nAttempts=0
        while i<nControls:
            newPoint = random.sample(self.points, 1)[0]
            criteriaFitters = 0
            suitable = True
            for c in controls:
                distance = osrm.getRoutes(newPoint,c).getDistances()[0]
                if distance < minDistance:
                    suitable = False
                    break
                elif distance < maxDistance:
                    criteriaFitters +=1
            if suitable and i>0 and criteriaFitters < math.floor(1+3*math.exp(-nAttempts/100)):
                suitable = False
            nAttempts+=1
            if suitable:
                print ("Control "+str(i)+" "+str(controlOrder[i])+" "+str(nAttempts))
                controls.append(Control(newPoint,controlOrder[i]))
                i+=1
                nAttempts=0
        course = Course(self.startFinish)
        course.addControls(controls)
        return course

class Parameters:
    '''Holds the limits of the course region, and the start/finish.'''
    def __init__(self, parameterFilename):
        tree = ET.parse(parameterFilename)
        ns = {"kml": KML_NAMESPACE}
        for placemark in tree.findall(".//kml:Placemark", ns):
            name = placemark.find("kml:name", ns).text
            if name.lower() == "limit":
                ring = [tuple([float(coord) for coord in coordinate.split(",")]) for coordinate in placemark.find(".//kml:coordinates", ns).text.split()]
                self.limit = geometry.Polygon(ring)
            elif name.lower() == "start":
                startCoords = [float(coord) for coord in placemark.find(".//kml:coordinates", ns).text.split(",")]
                self.start = Point(startCoords[1],startCoords[0])

    def restrict(self, pointSet):
        '''Remove all points outside the desired region.'''
        return [point for point in pointSet if self.limit.contains(geometry.Point((point.lon,point.lat)))]

if __name__ == "__main__":
    f=open(r"C:\Users\TC\Orienteering\Auto\lightLocations.csv")
    points = [tuple(line.rstrip().split(",")) for line in f]
    pointSet = [Point.fromBNG(point[0],point[1]) for point in points]
    params = Parameters(r"C:\Users\TC\Orienteering\Auto\limit_test.kml")
    restrictedPoints = params.restrict(pointSet)
    candidates = Candidates(restrictedPoints,params.start)
    course = candidates.courseByAddition(30, 400, 700)
    #course = candidates.courseByMinDistance(30,300)
    #course = candidates.courseByTwistiness(30)

    '''i=1
    controls = random.sample(restrictedPoints,1)
    while i<30:
        newPoint = random.choice(restrictedPoints)
        suitable = True
        for c in controls:
            r = requests.get(
                "http://192.168.1.85:5000/route/v1/driving/" + newPoint[1] + "," + newPoint[0] + ";" +
                c[1] + "," + c[0], params={"steps": "true", "alternatives": "true"})
            rj = r.json()
            if rj['routes'][0]['distance']<250:
                suitable = False
                break
        if suitable:
            i+=1
            controls.append(newPoint)
    i=1
    for c in controls:
        course.addControls([Control(Point(c[0],c[1]),i)])
        i+=1
    i=1
    while True:
        pointPair = random.sample(restrictedPoints,2)
        r = requests.get(
            "http://192.168.1.85:5000/route/v1/driving/" + pointPair[0][1] + "," + pointPair[0][0] + ";" + pointPair[1][
                1] + "," + pointPair[1][0], params={"steps": "true", "alternatives": "true"})
        rj = r.json()
        routes = rj['routes']
        if len(routes)>1:
            if routes[1]['distance']/routes[0]['distance']>1.1 and 1.0*len(routes[1]['legs'][0]['steps'])/len(routes[0]['legs'][0]['steps'])<0.6:
                print(i)
                course.addControls([Control(Point(pointPair[0][0], pointPair[0][1]), i)])
                i+=1
                course.addControls([Control(Point(pointPair[1][0], pointPair[1][1]), i)])
                i+=1
        if i>29:
            break
    for newPoint in random.sample(restrictedPoints,30):
        course.addControls([Control(Point(newPoint[0],newPoint[1]),i)])
        i+=1'''
    osrm = Osrm(OSRM_SERVER)
    print("Roundtrip ",osrm.getRoundtripLength(course))
    courseRep = CourseRepresentation(course)
    courseRep.makeKmlFile(r"C:\Users\TC\Orienteering\Auto\test.kml")
    courseRep.makeKmzFile(r"C:\Users\TC\Orienteering\Auto\test.kmz","Title")
    #courseRep.makePpenFile()