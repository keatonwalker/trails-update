"""Tool to produce information about trail geometery."""
import arcpy
import os
import matplotlib.pyplot as plt
from time import strftime


# Give FieldMap a new method to aviod the annoying output name pattern
def setOutputFieldName(self, name):
    """set the output field name of a fieldmap."""
    tempName = self.outputField
    tempName.name = name
    tempName.aliasName = name
    self.outputField = tempName
arcpy.FieldMap.outputFieldSettyBetty = setOutputFieldName

uniqueRunNum = strftime("%Y%m%d_%H%M%S")
outputWorkspace = r'C:\gis_working\fs trails\outputs.gdb'
outputFields = (('SHAPE@', 'geometery'), ('angle', 'FLOAT'),
                ('percent_len', 'FLOAT'), ('dist1', 'FLOAT'),
                ('heading1', 'FLOAT'), ('dist2', 'FLOAT'),
                ('heading2', 'FLOAT'))


class Feature (object):
    """Store usefull feature class information."""

    def __init__(self, workspace, name, spatialRef=None):
        """constructor."""
        self.workspace = workspace
        self.name = name
        self.path = os.path.join(workspace, name)
        self.spatialReference = spatialRef
        if self.spatialReference is None:
            self.spatialReference = arcpy.Describe(self.path).spatialReference

    def createThisFeature(self, goeType, fieldList=[]):
        """Create a feature with instance fields."""
        pass

    @staticmethod
    def createFeature(workspace, name, spatialRef, geoType, fieldList=[]):
        """Create a feature class and retrun a feature object."""
        tempFeature = Feature(workspace, name, spatialRef)
        arcpy.CreateFeatureclass_management(workspace,
                                            name,
                                            geoType,
                                            spatial_reference=spatialRef)
        if len(fieldList) > 0:
            for field in fieldList:
                name = field[0]
                fieldType = field[1]
                if name == 'SHAPE@':
                    continue
                arcpy.AddField_management(tempFeature.path,
                                          name,
                                          fieldType)

        return tempFeature


def plotAngleDist(dists, angles):
    """plot angle and distance values."""
    plt.plot(dists, angles)
    plt.show()


def getAngleDistStats(trailsFeature):
    """Get angle and distance stats for a road feature."""
    angles = []
    vertexDists = []
    outputVertexRows = []
    spatialRef = trailsFeature.spatialReference
    outputGdb = r"C:\gis_working\fs trails\Temp.gdb"
    outputVertexName = 'vertex_' + uniqueRunNum
    anglePoints = Feature.createFeature(outputGdb, outputVertexName,
                                        spatialRef, 'POINT', outputFields)
    outputVertexFeature = anglePoints.path

    with arcpy.da.SearchCursor(trailsFeature.path,
                               ['SHAPE@']) as lineCursor:
        for lineRow in lineCursor:
            line = lineRow[0]
            points = line.getPart(0)
            i = 0
            while i < line.pointCount - 2:
                end1 = arcpy.PointGeometry(points[i], spatialRef)
                vertex = arcpy.PointGeometry(points[i+1], spatialRef)
                end2 = arcpy.PointGeometry(points[i+2], spatialRef)
                heading1, dist1 = vertex.angleAndDistanceTo(end1)
                heading2, dist2 = vertex.angleAndDistanceTo(end2)
                if heading1 < 0:
                    heading1 = heading1 + 360
                if heading2 < 0:
                    heading2 = heading2 + 360
                angle = abs(heading1 - heading2)
                if angle > 180:
                    angle = 360 - angle
                angles.append(angle)
                vertexDist = line.measureOnLine(vertex, True)
                vertexDists.append(vertexDist)

                outputVertexRows.append((vertex, angle, vertexDist,
                                         dist1, heading1, dist2, heading2))
                i += 1

        vertexCursor = arcpy.da.InsertCursor(outputVertexFeature,
                                             [x[0] for x in outputFields])
        for row in outputVertexRows:
            vertexCursor.insertRow(row)
        del vertexCursor

        return (vertexDists, angles)


def getRenameFieldMap(featurePath, currentName, newName):
    """Create a field map that does a basic field rename."""
    tempMap = arcpy.FieldMap()
    tempMap.addInputField(featurePath, currentName)

    tempName = tempMap.outputField
    tempName.name = newName
    tempName.aliasName = newName
    tempMap.outputField = tempName
    # tempMap.outputFieldSettyBetty(newName)

    return tempMap


def featurePointCompare(srcLines, newLines):
    """Convert lines to points and compare position and other stats."""
    # Add lineLength field that will be present in converted point
    # Add lineId field to keep track of OBJECTID and aviod ESRI renaming
    lengthField = 'lineLength'
    lineIdField = 'lineId'
    newLengthField = 'nlineLength'  # Unique field names needed for field map.
    newLineIdField = 'nlineId'
    arcpy.AddField_management(srcLines.path, lengthField, 'DOUBLE')
    arcpy.CalculateField_management(srcLines.path,
                                    lengthField,
                                    "!shape.length@METERS!",
                                    'PYTHON_9.3')
    arcpy.AddField_management(srcLines.path, lineIdField, 'LONG')
    arcpy.CalculateField_management(srcLines.path,
                                    lineIdField,
                                    "!OBJECTID!",
                                    'PYTHON_9.3')

    arcpy.AddField_management(newLines.path, newLengthField, 'DOUBLE')
    arcpy.CalculateField_management(newLines.path,
                                    newLengthField,
                                    "!shape.length@METERS!",
                                    'PYTHON_9.3')
    arcpy.AddField_management(newLines.path, newLineIdField, 'LONG')
    arcpy.CalculateField_management(newLines.path,
                                    newLineIdField,
                                    "!OBJECTID!",
                                    'PYTHON_9.3')
    # Convert lines to a centroid point that is on the line.
    srcPoints = Feature(outputWorkspace,
                        'src_FeatToPt_' + uniqueRunNum,
                        srcLines.spatialReference)
    arcpy.FeatureToPoint_management(srcLines.path,
                                    srcPoints.path,
                                    'INSIDE')
    newPoints = Feature(outputWorkspace,
                        'new_FeatToPt_' + uniqueRunNum,
                        newLines.spatialReference)
    arcpy.FeatureToPoint_management(newLines.path,
                                    newPoints.path,
                                    'INSIDE')
    # Create field mappings
    pntJoinFMs = arcpy.FieldMappings()
    # Use field mapping to get line length difference
    lengthFM = arcpy.FieldMap()
    lengthFM.addInputField(srcPoints.path, lengthField)
    lengthFM.addInputField(newPoints.path, newLengthField)
    lengthFM.mergeRule = 'Range'
    lengthFM.outputFieldSettyBetty('lengthDiff')
    # Perform some field renaming
    srcIdFM = getRenameFieldMap(srcPoints.path, lineIdField, 'srcId')
    srcLengthFM = getRenameFieldMap(srcPoints.path, lengthField, 'srcLength')
    newIdFM = getRenameFieldMap(newPoints.path, newLineIdField, 'newId')
    newLengthFM = getRenameFieldMap(newPoints.path, newLengthField, 'newLength')
    # Add order determines table field order
    pntJoinFMs.addFieldMap(srcIdFM)
    pntJoinFMs.addFieldMap(newIdFM)
    pntJoinFMs.addFieldMap(srcLengthFM)
    pntJoinFMs.addFieldMap(newLengthFM)
    pntJoinFMs.addFieldMap(lengthFM)
    # Join the src (SGID based) points to the new (NFS based) points
    joinPoints = Feature(outputWorkspace,
                         'centriodPntJoin_' + uniqueRunNum,
                         srcPoints.spatialReference)
    arcpy.SpatialJoin_analysis(newPoints.path,
                               srcPoints.path,
                               joinPoints.path,
                               field_mapping=pntJoinFMs,
                               match_option='CLOSEST',
                               distance_field_name='joinDistance')

    lengthPercentageField = 'lenDiffPercent'
    arcpy.AddField_management(joinPoints.path, lengthPercentageField, 'FLOAT')
    arcpy.CalculateField_management(joinPoints.path,
                                    lengthPercentageField,
                                    "!lengthDiff!/!newLength! * 100",
                                    'PYTHON_9.3')

    joinDistField = 'joinDistPercent'
    arcpy.AddField_management(joinPoints.path, joinDistField, 'FLOAT')
    arcpy.CalculateField_management(joinPoints.path,
                                    joinDistField,
                                    "!joinDistance!/!newLength! * 100",
                                    'PYTHON_9.3')


if __name__ == '__main__':
    print 'hello'
    dataGdb = r'C:\gis_working\fs trails\Temp.gdb'
    trailsFeature = Feature(dataGdb,
                            'One_NFS_trail')

    # x, y = getAngleDistStats(trailsFeature)
    # plotAngleDist(x, y)
    srcLines = Feature(dataGdb, 'SGID_Full')
    newLines = Feature(dataGdb, 'TrailNFS_full_Project')
    featurePointCompare(srcLines, newLines)
