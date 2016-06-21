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
# outputWorkspace = r'C:\gis_working\fs trails\OutputsNfs.gdb'
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


def plotAngleDist(dists, angles, smaAngles, emaAngles):
    """plot angle and distance values."""
    plt.figure(1)
    plt.subplot(311)
    plt.plot(dists, angles, 'b')
    # plt.plot(dists, angles, 'bs')
    plt.subplot(312)
    plt.plot(dists, smaAngles, 'g')
    # plt.plot(dists, smaAngles, 'g^')
    plt.subplot(313)
    plt.plot(dists, emaAngles, 'y')
    plt.show()


def calcMovingAverage(currentAverage, newValue, lastValue, n):
    """Simple moving average helper."""
    return currentAverage + (newValue / float(n)) - (lastValue / float(n))


def calcExponetialMovingAve(currentAverage, newValue, alpha):
    """EMA helper."""
    return alpha * newValue + (1 - alpha) * currentAverage


def getAngleDistStats(trailsFeature):
    """Get angle and distance stats for a road feature."""
    smaAngles = []
    emaAngles = []
    emaAlpha = 0.20
    smaWindow = 20
    angles = []
    vertexDists = []
    outputVertexRows = []
    spatialRef = trailsFeature.spatialReference
    outputGdb = outputWorkspace
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
                # Calculate the simple moving average of angles
                if i + 1 <= smaWindow:
                    newSma = sum(angles) / float(len(angles))
                    smaAngles.append(newSma)
                else:
                    newSma = calcMovingAverage(smaAngles[i - 1], angle, smaAngles[i - smaWindow], smaWindow)
                    smaAngles.append(newSma)
                # Calculate exponential moving average of angles
                if i == 0:
                    emaAngles.append(angle)
                else:
                    newEma = calcExponetialMovingAve(emaAngles[i - 1], angle, emaAlpha)
                    emaAngles.append(newEma)

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

        return (vertexDists, angles, smaAngles, emaAngles)


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

    resultField = 'resultCategory'
    arcpy.AddField_management(joinPoints.path, resultField, 'TEXT', field_length=30)

    sameLayer = 'same'
    sameSelection = '{} < 10 and {} < 3'.format(lengthPercentageField, joinDistField)
    arcpy.MakeFeatureLayer_management(joinPoints.path, sameLayer, sameSelection)
    arcpy.CalculateField_management(sameLayer,
                                    resultField,
                                    '"same"',
                                    'PYTHON_9.3')

    additionLayer = 'additions'
    distFromSrcLines = 100
    additionSelection = '{} > 20'.format(joinDistField)
    arcpy.MakeFeatureLayer_management(joinPoints.path, additionLayer, additionSelection)
    arcpy.SelectLayerByLocation_management(additionLayer,
                                           'WITHIN_A_DISTANCE',
                                           srcLines.path,
                                           distFromSrcLines,
                                           invert_spatial_relationship=True)
    arcpy.CalculateField_management(additionLayer,
                                    resultField,
                                    '"addition"',
                                    'PYTHON_9.3')
    return joinPoints


def getExtendPoints(lineEndPoint, targetPoint, referenceLine, targetLine):
    """Extend referenceLine to the targetLine."""
    touchLine = arcpy.Polyline(arcpy.Array([lineEndPoint.centroid, targetPoint.centroid]),
                               referenceLine.spatialReference)
    extendPoints = []
    pnt33 = touchLine.positionAlongLine(0.33, True)
    pntRef33 = referenceLine.snapToLine(pnt33)
    pnt66 = touchLine.positionAlongLine(0.66, True)
    pntRef66 = referenceLine.snapToLine(pnt66)
    pnt90 = touchLine.positionAlongLine(0.90, True)
    pntRef90 = referenceLine.snapToLine(pnt90)
    newTargetPnt = targetLine.snapToLine(pntRef90)
    extendPoints = [lineEndPoint, pntRef33, pntRef66, pntRef90, newTargetPnt]
    for p in extendPoints:
        if lineEndPoint.distanceTo(p) > 100:  # remove points if they moved far away
            return [lineEndPoint, newTargetPnt]
    return extendPoints


def addNewLines(comparePoints, newLines, srcLines):
    """Process addition lines and connect ones that are within a distance of source lines."""
    newId = 'newId'
    resultField = 'resultCategory'
    additionCategory = 'addition'
    disconectedCategory = 'notouch_addition'
    conectedCategory = 'touch_addition'
    touchingDist = 50
    # Join result categories to the new lines
    arcpy.JoinField_management(newLines.path, 'OBJECTID', comparePoints.path, newId, [resultField])
    additionLayer = 'addLines'
    additionSelection = "{} = '{}'".format(resultField, additionCategory)
    arcpy.MakeFeatureLayer_management(newLines.path, additionLayer, additionSelection)
    arcpy.SelectLayerByLocation_management(additionLayer,
                                           'WITHIN_A_DISTANCE',
                                           srcLines.path,
                                           touchingDist,
                                           invert_spatial_relationship=True)
    arcpy.CalculateField_management(additionLayer,
                                    resultField,
                                    '"{}"'.format(disconectedCategory),
                                    'PYTHON_9.3')

    arcpy.SelectLayerByLocation_management(additionLayer,
                                           'WITHIN_A_DISTANCE',
                                           srcLines.path,
                                           touchingDist,
                                           invert_spatial_relationship=False)
    arcpy.CalculateField_management(additionLayer,
                                    resultField,
                                    '"{}"'.format(conectedCategory),
                                    'PYTHON_9.3')
    # Erase newLine parts that are within touchingDist
    srcLayer = 'srcLines'
    arcpy.MakeFeatureLayer_management(srcLines.path, srcLayer)
    arcpy.SelectLayerByLocation_management(srcLayer,
                                           'WITHIN_A_DISTANCE',
                                           additionLayer,  # newlines with touch_addition still selected
                                           touchingDist,
                                           invert_spatial_relationship=False)
    eraseBuffer = Feature(outputWorkspace, 'srcEraseBuff_' + uniqueRunNum, srcLines.spatialReference)
    arcpy.Buffer_analysis(srcLayer, eraseBuffer.path, touchingDist)
    eraseBufferLayer = 'eraseBuffer'
    arcpy.MakeFeatureLayer_management(eraseBuffer.path, eraseBufferLayer)
    touchLayer = 'touchLines'
    touchSelection = "{} = '{}'".format(resultField, 'touch_addition')
    arcpy.MakeFeatureLayer_management(newLines.path, touchLayer, touchSelection)
    erasedLines = Feature(outputWorkspace, 'erasedLines_' + uniqueRunNum, newLines.spatialReference)
    arcpy.Erase_analysis(touchLayer, eraseBuffer.path, erasedLines.path)
    erasedSingleLines = Feature(outputWorkspace, 'erasedSinglePart_' + uniqueRunNum, newLines.spatialReference)
    arcpy.MultipartToSinglepart_management(erasedLines.path, erasedSingleLines.path)
    # Update ereased lines to touch srcLines
    srcTouchLayer = 'srcTouch'
    arcpy.MakeFeatureLayer_management(srcLines.path, srcTouchLayer)
    arcpy.SelectLayerByLocation_management(srcTouchLayer,
                                           'WITHIN_A_DISTANCE',
                                           erasedSingleLines.path,
                                           touchingDist)
    srcCursor = arcpy.da.SearchCursor(srcTouchLayer, 'SHAPE@')
    srcTouchLines = [row[0] for row in srcCursor]
    del srcCursor
    i = 1
    with arcpy.da.UpdateCursor(erasedSingleLines.path, 'SHAPE@') as eraseCursor:
        for row in eraseCursor:
            updateLine = row[0]
            start = arcpy.PointGeometry(row[0].firstPoint, newLines.spatialReference)
            end = arcpy.PointGeometry(row[0].lastPoint, newLines.spatialReference)
            startPoints = []
            endPoints = []
            # Get the original geometry for the erase line
            arcpy.SelectLayerByLocation_management(touchLayer,
                                                   'INTERSECT',
                                                   start,
                                                   1)
            refLine = arcpy.CopyFeatures_management(touchLayer, arcpy.Geometry())
            for sl in srcTouchLines:
                if start is None and end is None:
                    break
                if start is not None and start.distanceTo(sl) <= touchingDist + 2:  # + 1 for tolerance
                    i += 1
                    startPoints.extend(getExtendPoints(start, sl.snapToLine(start), refLine[0], sl))
                    startPoints.reverse()
                    # startPointsDebug.append(sl.snapToLine(start))
                    start = None
                if end is not None and end.distanceTo(sl) <= touchingDist + 2:
                    endPoints.extend(getExtendPoints(end, sl.snapToLine(end), refLine[0], sl))
                    # startPointsDebug.append(sl.snapToLine(end))
                    i += 1
                    end = None
            if len(startPoints) > 0 or len(endPoints) > 0:
                updatePoints = updateLine.getPart(0)
                startPoints = [p.centroid for p in startPoints]
                endPoints = [p.centroid for p in endPoints]
                startPoints.extend(list(updatePoints))
                startPoints.extend(endPoints)
                uLine = arcpy.Polyline(arcpy.Array(startPoints), updateLine.spatialReference)
                row[0] = uLine
                eraseCursor.updateRow(row)

    print i
    # startPointsDebug = [arcpy.PointGeometry(p, newLines.spatialReference) for p in startPointsDebug]
    # arcpy.CopyFeatures_management(startPointsDebug, os.path.join(outputWorkspace, 'tempStart_' + uniqueRunNum))

    disconnectedLayer = 'disconectedAdd'
    disSelection = "{} = '{}'".format(resultField, disconectedCategory)
    arcpy.MakeFeatureLayer_management(newLines.path, disconnectedLayer, disSelection)
    arcpy.Append_management(disconnectedLayer, erasedSingleLines.path, 'NO_TEST')

    arcpy.Delete_management(additionLayer)
    arcpy.Delete_management(disconnectedLayer)
    arcpy.Delete_management(srcLayer)
    arcpy.Delete_management(eraseBufferLayer)
    arcpy.Delete_management(touchLayer)

    return erasedSingleLines


def mergeAdditionLines(srcLines, addLines, fieldMap, updateSource):
    """Merge and addLines and srcLines and map usefull fields."""
    dataSourceField = 'UpdateSource'
    arcpy.AddField_management(srcLines.path, dataSourceField, 'TEXT')
    arcpy.AddField_management(addLines.path, dataSourceField, 'TEXT')
    arcpy.CalculateField_management(addLines.path,
                                    dataSourceField,
                                    "'{} !{}!'".format(updateSource, 'OBJECTID'),
                                    'PYTHON_9.3')
    updateSrcFM = getRenameFieldMap(addLines.path, dataSourceField, dataSourceField)
    fieldMap.addFieldMap(updateSrcFM)

    arcpy.Append_management(addLines.path, srcLines.path, 'NO_TEST', fieldMap)


def getCacheFieldMap(addLines):
    """Get field map for Cache."""
    srcNameField = 'PrimaryName'
    addNameField = 'primarynam'
    # Create field mappings
    additionFMs = arcpy.FieldMappings()
    # Perform some field renaming
    nameFM = getRenameFieldMap(addLines.path, addNameField, srcNameField)
    trailIdFM = getRenameFieldMap(addLines.path, 'trialid', 'TrlID')
    descFM = getRenameFieldMap(addLines.path, 'descriptio', 'Description')
    desgUsesFM = getRenameFieldMap(addLines.path, 'designated', 'DesignatedUses')
    surfTypeFM = getRenameFieldMap(addLines.path, 'surfacetyp', 'surfacetype')
    motoFM = getRenameFieldMap(addLines.path, 'motorizedp', 'MotorizedProhibited')
    maintainFM = getRenameFieldMap(addLines.path, 'primarymai', 'PrimaryMaint')
    adaFM = getRenameFieldMap(addLines.path, 'adaaccessi', 'ADAAccessible')

    additionFMs.addFieldMap(trailIdFM)
    additionFMs.addFieldMap(descFM)
    additionFMs.addFieldMap(desgUsesFM)
    additionFMs.addFieldMap(nameFM)
    additionFMs.addFieldMap(surfTypeFM)
    additionFMs.addFieldMap(motoFM)
    additionFMs.addFieldMap(maintainFM)
    additionFMs.addFieldMap(adaFM)

    return additionFMs


def getUtahFieldMap(addLines):
    """Get field map for Utah."""
    srcNameField = 'PrimaryName'
    addNameField = 'Name'
    # Create field mappings
    additionFMs = arcpy.FieldMappings()
    # Perform some field renaming
    nameFM = getRenameFieldMap(addLines.path, addNameField, srcNameField)
    sourceFM = getRenameFieldMap(addLines.path, 'SourceCity', 'DataSource')
    surfTypeFM = getRenameFieldMap(addLines.path, 'Surface', 'SurfaceType')
    commentFM = getRenameFieldMap(addLines.path, 'Category1', 'Comments')

    additionFMs.addFieldMap(nameFM)
    additionFMs.addFieldMap(sourceFM)
    additionFMs.addFieldMap(surfTypeFM)
    additionFMs.addFieldMap(commentFM)

    return additionFMs


def getWasatchFieldMap(addLines):
    """Get field map for Wasatch."""
    srcNameField = 'PrimaryName'
    addNameField = 'NAME'
    # Create field mappings
    additionFMs = arcpy.FieldMappings()
    # Perform some field renaming
    nameFM = getRenameFieldMap(addLines.path, addNameField, srcNameField)
    surfTypeFM = getRenameFieldMap(addLines.path, 'Surface', 'SurfaceType')
    desgUsesFM = getRenameFieldMap(addLines.path, 'Use_', 'DesignatedUses')
    maintainFM = getRenameFieldMap(addLines.path, 'Jurisdicti', 'PrimaryMaint')
    statusFM = getRenameFieldMap(addLines.path, 'Status', 'Status')

    additionFMs.addFieldMap(nameFM)
    additionFMs.addFieldMap(surfTypeFM)
    additionFMs.addFieldMap(desgUsesFM)
    additionFMs.addFieldMap(maintainFM)
    additionFMs.addFieldMap(statusFM)

    return additionFMs


def getNfsFieldMap(addLines):
    """Get field map for Wasatch."""
    srcNameField = 'PrimaryName'
    addNameField = 'TRAIL_NAME'
    # Create field mappings
    additionFMs = arcpy.FieldMappings()
    # Perform some field renaming
    nameFM = getRenameFieldMap(addLines.path, addNameField, srcNameField)
    surfTypeFM = getRenameFieldMap(addLines.path, 'TRAIL_TYPE', 'SurfaceType')
    trailIdFM = getRenameFieldMap(addLines.path, 'TRAIL_NO', 'TrailID')
    sourceFM = getRenameFieldMap(addLines.path, 'ATTRIBUTESUBSET', 'DataSource')

    additionFMs.addFieldMap(nameFM)
    additionFMs.addFieldMap(surfTypeFM)
    additionFMs.addFieldMap(trailIdFM)
    additionFMs.addFieldMap(sourceFM)

    return additionFMs


if __name__ == '__main__':
    print uniqueRunNum
    global outputWorkspace
    outputWorkspace = r'C:\gis_working\fs trails\OutputsNfs.gdb'
    dataGdb = r'C:\gis_working\fs trails\sourceData.gdb'
    # trailsFeature = Feature(dataGdb, 'One_NFS_trail')

    # x, y, sma, ema = getAngleDistStats(trailsFeature)
    # plotAngleDist(x, y, sma, ema)

    baseSrcLines = Feature(dataGdb, 'SGID_Cache_Utah_Was')
    srcLines = Feature(outputWorkspace, baseSrcLines.name + '_' + uniqueRunNum, baseSrcLines.spatialReference)
    arcpy.CopyFeatures_management(baseSrcLines.path, srcLines.path)
    baseNewLines = Feature(dataGdb, 'NFS_full_base')
    newLines = Feature(outputWorkspace, baseNewLines.name + '_' + uniqueRunNum, baseNewLines.spatialReference)
    arcpy.CopyFeatures_management(baseNewLines.path, newLines.path)

    comparePoints = featurePointCompare(srcLines, newLines)
    additionLines = addNewLines(comparePoints, newLines, srcLines)
    mergeAdditionLines(srcLines, additionLines, getNfsFieldMap(additionLines), 'NFS')
