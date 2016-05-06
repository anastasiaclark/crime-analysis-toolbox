# -----------------------------------------------------------------------------
# Copyright 2016 Esri
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# -----------------------------------------------------------------------------

# ==================================================
# calculate_prediction_zones.py BETA
# --------------------------------------------------
# requirments: ArcMap/ArcCatalog 10.3.1+ with Spatial Analyst,
#              ArcGIS Pro 1.2+ with Spatial Analyst
#              ArcGIS Standard or Advanced license required
#              Python 2.7 or 3.4
# author: ArcGIS Solutions
# contact: ArcGISTeamLocalGov@esri.com
# company: Esri
# ==================================================
# description: Calculate a probability surface showing the liklihood of the
#              occurance of a repeat or near repeat incident
# ==================================================
# history:
# v1 03/25/2016
# ==================================================

import arcpy

from datetime import datetime as dt
from datetime import date as dy
from datetime import timedelta as td
from os import path
import os
import math

from arcrest.security import AGOLTokenSecurityHandler
from arcresthelper import securityhandlerhelper
from arcresthelper import common
from arcrest.agol import FeatureLayer

# Enable overwriting datasets
arcpy.env.overwriteOutput = True

# Status fields for output polygons (created if non-existant)
cur_status_field = 'MOSTRECENT'
cur_date_field = 'CREATEDATE'
proc_date_field = 'STARTDATE'
temporal_band_field = 'TIMEBAND'
spatial_band_field = 'SPACEBAND'
risk_range_field = "RISKRANGE"

# Get current date & time
today = dy.today()
todaytime = dt.today()


def expand_extents(data, stretch):
    """Expand the extents of a dataset by a set distance in all directions by
        adding and removing features from the dataset"""

    # Calculate new dataset extents from current values
    d = arcpy.Describe(data)

    xmin = d.extent.XMin - stretch
    ymax = d.extent.YMax + stretch

    new_features = [(xmin, d.extent.YMin - stretch),
                    (d.extent.XMax + stretch, ymax)]

    # Add two features to the dataset at the new xmin, ymin and xmax, ymax
    with arcpy.da.InsertCursor(data, ['Shape@X', 'Shape@Y']) as cursor:
        for feat in new_features:
            cursor.insertRow(feat)

    # Delete the two new features based on their coordinate values
    with arcpy.da.UpdateCursor(data, ['Shape@X', 'Shape@Y']) as rows:
        for row in rows:
            if row[0] == xmin or row[1] == ymax:
                rows.deleteRow()

# End of expand_extents function


def connect_to_layer(username, password, server_url, service_url):
    """Connect to ArcGIS Online or Portal for ArcGIS layer"""
    proxy_port = None
    proxy_url = None

    si = {}
    si['security_type'] = 'Portal'  # LDAP, NTLM, OAuth, Portal, PKI
    si['username'] = username
    si['password'] = password
    si['org_url'] = server_url
    si['proxy_url'] = proxy_url
    si['proxy_port'] = proxy_port
    si['referer_url'] = None
    si['token_url'] = None
    si['certificatefile'] = None
    si['keyfile'] = None
    si['client_id'] = None
    si['secret_id'] = None

    shh = securityhandlerhelper.securityhandlerhelper(securityinfo=si)
    if not shh.valid:
        raise Exception(shh.message)

    fl = FeatureLayer(
        url=service_url,
        securityHandler=shh.securityhandler,
        proxy_port=proxy_port,
        proxy_url=proxy_url,
        initialize=True)

    return fl

# End of connect_to_layer function


def calculate_risk_surface(lyr, age, dist, halflife, halfdist):
    """Create raster of risk extent and intensity based on incident age and
       spatial reach of influence"""

    # Build float distance raster for incident
    dist_raster = arcpy.sa.EucDistance(lyr,
                                       dist)

    # Apply distance & temporal decay - math.log() is ln()
    max_temporal_risk = 1.0

    k_temporal = math.log(0.5)/-halflife
    k_spatial = math.log(0.5)/-halfdist

    dist_raster = arcpy.sa.Exp(dist_raster * -k_spatial)

    inc_raster = max_temporal_risk * math.exp(-age * k_temporal) * dist_raster

    # Set Null values to 0 to allow for raster math
    null_locations = arcpy.sa.IsNull(inc_raster)
    inc_raster = arcpy.sa.Con(null_locations, 0,
                              inc_raster, where_clause="Value = 1")
    return inc_raster

# End of calculate_risk_surface function


def calculate_max_risk(cumulative_raster, new_raster):
    """Creates a raster using the maximum values of two rasters"""
    # Determine cells where sum raster has => value
    sum_vals = arcpy.sa.GreaterThanEqual(cumulative_raster, new_raster)
    sum_vals = arcpy.sa.Times(sum_vals, cumulative_raster)

    # Determine cells where inc raster has greater value
    inc_vals = arcpy.sa.GreaterThan(new_raster, cumulative_raster)
    inc_vals = arcpy.sa.Times(inc_vals, new_raster)

    # Sum greatest value rasters
    cumulative_raster += inc_vals

    return cumulative_raster

# End of calculate_max_risk function


def add_status_fields_to_lyr(lyr):
    """Adds a text and/or date field to a fc or lyr for tracking status"""
    fields = [f.name for f in arcpy.ListFields(lyr)]

    if cur_status_field not in fields:
        arcpy.AddField_management(lyr, cur_status_field, 'TEXT', field_length=5)

    if cur_date_field not in fields:
        arcpy.AddField_management(lyr, cur_date_field, 'DATE')

    if risk_range_field not in fields:
        arcpy.AddField_management(lyr, risk_range_field, 'LONG')

    if proc_date_field not in fields:
        arcpy.AddField_management(lyr, proc_date_field, 'DATE')

    if temporal_band_field not in fields:
        arcpy.AddField_management(lyr, temporal_band_field, 'FLOAT')

    if spatial_band_field not in fields:
        arcpy.AddField_management(lyr, spatial_band_field, 'FLOAT')

# End of add_status_fields_to_lyr function


def add_status_field_to_service(fl):
    """Adds a text and/or date field to a hosted service for tracking status"""
    layer_fields = [f['name'] for f in fl.fields]
    fieldToAdd = {"fields": []}

    if cur_status_field not in layer_fields:
        fieldToAdd["fields"].append({
                "name": cur_status_field,
                "type": "esriFieldTypeString",
                "alias": cur_status_field,
                "sqlType": "sqlTypeOther",
                "length": 5,
                "nullable": True,
                "editable": True,
                "domain": None,
                "defaultValue": None})

    if cur_date_field not in layer_fields:
        fieldToAdd["fields"].append({
                "name": cur_date_field,
                "type": "esriFieldTypeDate",
                "alias": cur_date_field,
                "nullable": True,
                "editable": True,
                "domain": None,
                "defaultValue": None})

    if risk_range_field not in layer_fields:
        fieldToAdd["fields"].append({
                "name": risk_range_field,
                "type": "esriFieldTypeInteger",
                "alias": risk_range_field,
                "nullable": True,
                "editable": True,
                "domain": None,
                "defaultValue": None})

    if proc_date_field not in layer_fields:
        fieldToAdd["fields"].append({
                "name": proc_date_field,
                "type": "esriFieldTypeDate",
                "alias": proc_date_field,
                "nullable": True,
                "editable": True,
                "domain": None,
                "defaultValue": None})

    if temporal_band_field not in layer_fields:
        fieldToAdd["fields"].append({
                "name": temporal_band_field,
                "type": "esriFieldTypeDouble",
                "alias": temporal_band_field,
                "nullable": True,
                "editable": True,
                "domain": None,
                "defaultValue": None})

    if spatial_band_field not in layer_fields:
        fieldToAdd["fields"].append({
                "name": spatial_band_field,
                "type": "esriFieldTypeDouble",
                "alias": spatial_band_field,
                "nullable": True,
                "editable": True,
                "domain": None,
                "defaultValue": None})

    fl.administration.addToDefinition(fieldToAdd)

# End of add_status_field_to_service function


def convert_raster_to_zones(raster, bins, status_field, date_field, details):
    """Convert non-0 raster cell values to polygons using a
       set number of bins"""
    sliced = arcpy.sa.Slice(raster, int(bins))
    polys = arcpy.RasterToPolygon_conversion(sliced,
                                             path.join("in_memory",
                                                       "temp_polys"),
                                             "NO_SIMPLIFY")
    add_status_fields_to_lyr(polys)

    fields = [status_field, date_field, risk_range_field, "gridcode", proc_date_field, spatial_band_field, temporal_band_field]
    with arcpy.da.UpdateCursor(polys, fields) as rows:
        for row in rows:
            row[0] = 'True'
            row[1] = todaytime
            row[2] = row[3]
            row[4] = details[0]
            row[5] = details[1]
            row[6] = details[2]
            rows.updateRow(row)

    arcpy.DeleteField_management(polys, ["Id", "gridcode"])
    return polys

# End of convert_raster_to_zones function


def create_zone_fc(template, sr, out_path):
    """Create polygon feature class for prediction zone features"""
    poly_paths = out_path.split(os.sep)[:-1]
    poly_path = os.sep.join(poly_paths)
    poly_name = out_path.split(os.sep)[-1]

    arcpy.CreateFeatureclass_management(poly_path,
                                        poly_name,
                                        template=template)
    arcpy.DefineProjection_management(out_path, sr)

    return out_path

# End of create_zone_fc function


def main(in_features, date_field, init_date, spatial_band_size, spatial_half,
         temporal_band_size, temporal_half, probability_type, out_raster,
         out_polygon, slice_num, pub_polys='', pub_type='', username='',
         password='', server_url='', poly_url='', *args):

    """ Generates a raster and series of polygons based on that raster to
        illustrate the probability of incidents occuring at the current moment
        in time based on defined algorithms for the decay of spatial and
        temporal influence of previous incidents.

        in_features: Point feature class or shapefile showing the location of
                     incidents that have occured recently, and from which
                     predictions will be based. This feature class must have a
                     date field and all features must have date values.

        date_field: Field in in_features containing the date each incident
                    occurred. Values in this field are used to calculate the
                    decay of temporal influence between when the incident
                    occured and the current date.

        init_date: Initial processing date. All incidents with dates between
                   this date and init_date - temporal_band_size will be
                   included in the report

        spatial_band_size: Value in the units of in_features representing the
                           maximum reach of spatial influence of historical
                           incidents.

        spatial_half: Value between 0 and spatial_band_size representing the
                      distance at which the risk of a near repeat incident has
                      been halved. Can be taken from output of classification
                      tool.

        temporal_band_size: Value in days representing the maximum reach of
                            temporal influence of historical incidents.
                            Features in in_features where todays date minus the
                            incident date results in a number of days greater
                            than this value will not be considered when
                            creating the prediction zones.

        temporal_half: Value between 0 and temporal_band_size representing the
                      time at which the risk of a near repeat incident has
                      been halved. Can be taken from output of classification
                      tool.

        probability_type: 'CUMULATIVE' (default) creates a surface resulting
                          from summing the prediction risks from each incident;
                          'MAXIMUM' creates a surface representing the maximum
                          risk value from each incident.

        out_raster: Location for output incident prediction surface raster.
                    Raster name will have timestamp.

        out_polygon: Output polygon feature class based on classifying the
                        out_raster values into slice_num categories.
                        Polygon boundaries represent the bounds of the
                        prediction zones as defined by the raster slices.

        slice_num: Integer value representing the number of zones that will be
                   created from the prediction raster. Each zone will represent
                   a range of prediction risk values.

        pub_polys: booleen option for publishing the polygon features. Service
                   must exist previously. Service will be truncated and the
                   cumulative results from in_features will be appended

        init_date: initial processing date.
        pub_type: Choice of publication environments- NONE, ARCGIS_ONLINE,
                  ARCGIS_PORTAL, ARCGIS_SERVER

        username: administrative username for the service

        password: corresponding to the username

        server_url: organization url

        poly_url: URL to the rest endpoint of the polygon service layer
    """

    try:
        i = 0
        arcpy.SetProgressor("default")
        arcpy.SetProgressorLabel('Initializing...')

        # Check out spatial analyst extentsion
        if arcpy.CheckExtension("Spatial") == "Available":
            arcpy.CheckOutExtension("Spatial")
        else:
            raise Exception("Spatial Analyst license unavailable")

        now = dt.strftime(dt.now(), "%y%m%d%H%M%S")

        # Convert booleen values
        if not pub_polys == 'True':
            pub_polys = False

        # Get init_date value
        if init_date == 'TODAY':
            init_date = today
        elif init_date == 'YESTERDAY':
            init_date = today - td(days=1)
        else:
            try:
                init_date = dt.strptime(init_date, "%Y-%m-%d")
            except ValueError:
                raise Exception("Invalid date format. Initial Date must be in the format yyyy-mm-dd.")

        # test for input features
        in_count = arcpy.GetCount_management(in_features)
        if int(in_count.getOutput(0)) < 1:
            raise Exception("Insufficient features for processing")

        # Work in an in-memory copy of the dataset to avoid editing the original
        incident_fc = arcpy.FeatureClassToFeatureClass_conversion(in_features,
                                                                  "in_memory",
                                                                  'temp_incs')

        # Get OID field name
        oidname = arcpy.Describe(incident_fc).oidFieldName

        # Expand the extents of the dataset by the size of the spatial band
        #   rasters will represent the full extent of risk,
        #   not bound to extents of incidents
        expand_extents(incident_fc, float(spatial_band_size))

        # SelectLayerByAttributes tool requires feature layer
        incident_lyr = arcpy.MakeFeatureLayer_management(incident_fc)

        # Create in-memory summary raster with max extents
        d = arcpy.Describe(incident_fc)
        sr = d.spatialReference
        arcpy.env.extent = d.extent

        sum_raster = arcpy.sa.CreateConstantRaster(0, data_type='INTEGER',
                                                   extent=d.extent)

        # Calculate minimum bounds of accepted time frame
        date_min = init_date - td(days=int(temporal_band_size))

        # Create risk rasters for each incident within temporal reach of today
        sql = """{0} <= date'{1}' AND {0} > date'{2}'""".format(date_field,
                                                                init_date,
                                                                date_min)
        numrows = 0
        with arcpy.da.SearchCursor(incident_fc, "OID@", where_clause=sql) as rows:
            for row in rows:
                numrows += 1

        with arcpy.da.SearchCursor(incident_fc,
                                   ['OID@', date_field],
                                   where_clause=sql) as incidents:
            count = 0
            for incident in incidents:
                arcpy.SetProgressorLabel('Calculating influence of incident {} of {}...'.format(count+1, numrows))

                # Calculate age of incident
                try:
                    date_diff = init_date - incident[1].date()
                except TypeError:
                    date_diff = init_date.date() - incident[1].date()

                # Build float distance raster for incident
                sql = """{} = {}""".format(oidname, incident[0])
                arcpy.SelectLayerByAttribute_management(incident_lyr,
                                                        where_clause=sql)

                inc_raster = calculate_risk_surface(incident_lyr,
                                                    date_diff.days,
                                                    spatial_band_size,
                                                    float(temporal_half),
                                                    float(spatial_half))

                # Process cumulative risk
                if probability_type == 'CUMULATIVE':
                    sum_raster += inc_raster

                # Process maximum risk
                else:
                    sum_raster = calculate_max_risk(sum_raster, inc_raster)

                count += 1

        if not count:
            raise Exception('No incidents found between {} and {}'.format(date_min, init_date))
        else:
            arcpy.AddMessage("{} incidents found.".format(count))

        # Save final probability raster where values are > 0
        arcpy.SetProgressorLabel('Saving final raster...')

        sum_raster = arcpy.sa.SetNull(sum_raster, sum_raster, "Value <= 0")
        out_raster_name = ''.join([out_raster, os.sep, 'p', now])
        sum_raster.save(out_raster_name)
        arcpy.SetParameterAsText(18, out_raster_name)

        # Slice raster values into categories and convert to temp polys
        arcpy.SetProgressorLabel('Creating polygons...')

        temp_polys = convert_raster_to_zones(sum_raster, slice_num,
                                             cur_status_field, cur_date_field,
                                             [init_date, spatial_band_size,
                                             temporal_band_size])

        # Creat polygon fc if it doesn't exist
        if not arcpy.Exists(out_polygon):
            create_zone_fc(temp_polys, sr, out_polygon)

        # Create status fields if they don't exist
        add_status_fields_to_lyr(out_polygon)

        # Set status of all existing features to False
        sql = """{} <> 'False'""".format(cur_status_field)
        with arcpy.da.UpdateCursor(out_polygon,
                                   cur_status_field,
                                   where_clause=sql) as rows:
            for row in rows:
                row[0] = 'False'
                rows.updateRow(row)

        # Append temp poly features to output polygon fc
        arcpy.Append_management(temp_polys, out_polygon)
        arcpy.SetParameterAsText(17, out_polygon)

        # Update polygon services.
        # If pubtype = NONE or SERVER, no steps necessary

        if pub_type in ['ARCGIS_ONLINE', 'ARCGIS_PORTAL'] and pub_polys:
            arcpy.SetProgressorLabel('Updating polygon feature service...')

            # connect to incidents service
            try:
                fl = connect_to_layer(username, password, server_url, poly_url)
            except:
                raise Exception('Could not update service. Please verify '
                                'organization URL and service URL are '
                                'correct, and the provided username and '
                                'password have access to the service.')

            # Check service for status, creation, risk fields. add if necessary
            add_status_field_to_service(fl)

            # Update 'current' features in service to be 'past'
            field_info = [{'FieldName': cur_status_field,
                           'ValueToSet': 'False'}]

            out_fields = ['objectid']
            for fld in field_info:
                out_fields.append(fld['FieldName'])

            sql = """{} = 'True'""".format(cur_status_field)
            updateFeats = fl.query(where=sql,
                                   out_fields=','.join(out_fields))

            for feat in updateFeats:
                for fld in field_info:
                    feat.set_value(fld['FieldName'], fld['ValueToSet'])

            fl.updateFeature(features=updateFeats)

            # Add new 'current' features
            fl.addFeatures(temp_polys)

    except arcpy.ExecuteError:
        # Get the tool error messages
        msgs = arcpy.GetMessages()
        arcpy.AddError(msgs)
        print(msgs)

    except:

        # Return  error messages for use in script tool or Python Window
        arcpy.AddError(str(sys.exc_info()[1]))

        # Print Python error messages for use in Python / Python Window
        print("\n" + str(sys.exc_info()[1]) + "\n")

    finally:
        arcpy.SetProgressorLabel('Completed.')
        arcpy.CheckInExtension("Spatial")


if __name__ == '__main__':
    argv = [arcpy.GetParameterAsText(i)
            for i in range(arcpy.GetArgumentCount())]

    # Handle default values from results window
    i = 0
    while i < len(argv):
        if argv[i] == "#":
            argv[i] = ""
        i += 1
    argv = tuple(argv)

    main(*argv)
