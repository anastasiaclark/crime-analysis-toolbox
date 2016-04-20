# -----------------------------------------------------------------------------
# Copyright 2015 Esri
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
# incident_classification.py BETA
# --------------------------------------------------
# requirments: ArcMap/ArcCatalog 10.3.1+
#              ArcGIS Pro 1.2+
#              ArcGIS Advanced license required
#              Python 2.7 or 3.4
# author: ArcGIS Solutions
# contact: ArcGISTeamLocalGov@esri.com
# company: Esri
# ==================================================
# description: Classify incidents according to spatial and temporal proximity
#              to preceeding incidents. Generate a summary report of the
#              processed incidents
# ==================================================
# history:
# 04/19/2016 - AM - beta
# ==================================================

import arcpy
from datetime import datetime as dt
from datetime import timedelta as td
from os import path

# Added field names
spatial_band_field = 'SPACEBAND'
temporal_band_field = 'TIMEBAND'
incident_type_field = 'INCCLASS'
origin_feat_field = 'ORIGIN'
z_value_field = 'ZVALUE'
dist_orig_field = "DISTTOORIG"

units = {"Meter": "m",
         "Foot_US": "ft",
         "Foot": "ft",
         "150_Kilometers": "x 150km",
         "50_Kilometers": "x 50km"}


def reset_fields(fc):
    """Checks for required incident classification fields,
       and deletes/adds fields as necessary"""

    # Delete classification fields if they already exist in the dataset
    inc_fields = [f.name for f in arcpy.ListFields(fc)]

    delete_fields = []

    for field in[incident_type_field, z_value_field]:
        if field in inc_fields:
            delete_fields.append(field)

    if delete_fields:
        arcpy.DeleteField_management(fc, delete_fields)

    # Add field for incident classification
    arcpy.AddField_management(fc,
                              field_name=incident_type_field,
                              field_type='TEXT')

    # Add field for z value calculation
    arcpy.AddField_management(fc,
                              field_name=z_value_field,
                              field_type='LONG')


def calculate_band(value, bands):
    """Finds the first number in a list larger than a value"""
    for band in bands:
        if band > value:
            return band


def classify_incidents(in_features, date_field, report_location, repeatdist,
                       spatial_bands, temporal_bands, out_lines_dir,
                       out_lines_name, *args):
    """Updates an input feature class to classify features according to their
       proximity in space and time to previous incidents

       in_features: point feature class of incidents to classify. This dataset
                    will typically cover a large timespan. Must have a date
                    field.

       date_field: Field of type Date on the in_features dataset.
                   All features must have values in this field.

       report_location: Directory on disk where a summary report (csv) of the
                        processed incidents can be written.

       repeatdist: Distance in the units of in_features below which adjacent
                   incidents are considered repeats rather than near-repeats.
                   Default value is 0.

       spatial_bands: semi-colon separated list of distances in the unit of the
                      in_features. Features will be classified according to
                      the smallest value that exceeds their proximity in space
                      to the nearest preceeding incident that is also within
                      the maximum allowable temporal_band value.

       temporal_bands: semi-colon separated list of positive integers
                       representing the number of days between an originating
                       incident and a repeat or near repeat incident. Features
                       will be classified according to the smallest listed
                       value that exceeds their proximity in time to their
                       nearest spatial neighbour.

       out_lines_dir: The workspace where the line features will be stored

       out_lines_name: The name of the feature class that will be created to
                       hold the line features."""
    try:
        # Fix for potential issue with xlsx files as report locations
        if not path.isdir(report_location):
            report_location = path.dirname(report_location)

        # Build sorted lists of band values
        spatial_bands = [float(b) for b in spatial_bands.split(';')]
        temporal_bands = [float(b) for b in temporal_bands.split(';')]

        repeatdist = float(repeatdist)
        spatial_bands.append(repeatdist)

        spatial_bands = list(set(spatial_bands))
        temporal_bands = list(set(temporal_bands))

        spatial_bands.sort()
        temporal_bands.sort()

        arcpy.env.overwriteOutput = True

        # Report run time used for file names
        now = dt.strftime(dt.now(), "%Y-%m-%d_%H-%M-%S")
        now_nice = dt.strftime(dt.now(), "%Y-%m-%d %H:%M:%S")

        # Check for and delete existing fields necessary for classification
        reset_fields(in_features)

        # Get name of OID field
        oidname = arcpy.Describe(in_features).oidFieldName

        # Get sorted list of unique incident date values
        with arcpy.da.SearchCursor(in_features, date_field) as rows:
            date_vals = [row[0] for row in rows]

        date_vals = list(set(date_vals))
        date_vals.sort()

        # Range of incident dates
        min_date = date_vals[0]
        max_date = date_vals[-1]

        # Keep track of origins and nrs
        oids = []
        nrids = []
        rids = []

        # Connecting line segments and table rows
        new_lines = []
        new_rows = []

        # Build empty dictionary to hold spatial and temporal band tallies
        band_counts = {}
        for sband in spatial_bands:
            band_counts[sband] = {}
            for tband in temporal_bands:
                band_counts[sband][tband] = 0

        # Value lists for half life calculations
        all_distances = {}
        for sband in spatial_bands:
            all_distances[sband] = []

        all_lives = {}
        for tband in temporal_bands:
            all_lives[tband] = []

        found_connections = []

        # Build table of all records within the max spatial band of anther feature
        near_table = arcpy.GenerateNearTable_analysis(in_features, in_features, search_radius=temporal_bands[-1], closest='ALL', method='GEODESIC')

        # Identify and process relevent near features
        with arcpy.da.SearchCursor(near_table, field_names=['IN_FID', 'NEAR_FID', 'NEAR_DIST']) as nearrows:

            # Process each identified connection within the spatial bands
            for nearrow in nearrows:
                dist = nearrow[2]
                if not dist <= spatial_bands[-1]:
                    continue

                links= []

                # Find the two features that are part of the connection
                where_clause = """{} in ({},{})""".format(oidname, nearrow[0], nearrow[1])
                fields  = [oidname, date_field, incident_type_field, z_value_field, 'SHAPE@X','SHAPE@Y']
                with arcpy.da.UpdateCursor(in_features, field_names=fields, where_clause=where_clause) as cur_link:
                    for feat in cur_link:
                        # Calculate the z values of each incident in the pair
                        zval = feat[1] - min_date
                        feat[3] = zval.days
                        cur_link.updateRow(feat)
                        links.append([feat[0], feat[1], feat[4], feat[5], feat[3]])

                # Identify which feature is the oldest and id it as the source
                if links[0][1] > links[1][1]:
                    oid, odate, ox, oy, oz = links[1]
                    fid, fdate, fx, fy, fz = links[0]

                else:
                    oid, odate, ox, oy, oz = links[0]
                    fid, fdate, fx, fy, fz = links[1]

                # test for new connection
                if (oid, fid) in found_connections:
                    continue

                # Calculate the days between the two dates
                datediff = fdate - odate
                daydiff = datediff.days

                # only process rows within defined temporal bands
                if daydiff > temporal_bands[-1]:
                    continue

                # Identify the spatial bands that are covered by this relationship and create a connecting line feature
                link_found = False
                for sband in spatial_bands:
                    if dist <= sband and not link_found:
                        for tband in temporal_bands:
                            if daydiff <= tband and not link_found:
                                # count for band combo
                                band_counts[sband][tband] += 1

                                # track distances and lives for half measures
                                all_distances[sband].append(dist)
                                all_lives[tband].append(daydiff)
                                min_bands_found = True

                                # Repeats will only fall into first spatial band
                                if dist <= spatial_bands[0]:
                                    rids.append(fid)

                                link_found = True

                if link_found:
                    # track oid vs nrid
                    oids.append(oid)
                    nrids.append(fid)
                    found_connections.append((oid, fid))

                    # create connecting line from x, y, z values of two pts
                    end = arcpy.Point(X=fx, Y=fy, Z=fz)
                    start = arcpy.Point(X=ox, Y=oy, Z=oz)
                    vertices = arcpy.Array([start, end])
                    feature = arcpy.Polyline(vertices, None, True, False)
                    new_lines.append([fid, oid, dist, sband, daydiff, tband, feature])

        # Delete near table
        arcpy.Delete_management(near_table)

        # Create feature class for connecting lines
        sr = arcpy.Describe(in_features).spatialReference
        connectors = arcpy.CreateFeatureclass_management(out_lines_dir,
                                                         out_lines_name,
                                                         'POLYLINE',
                                                         has_z='ENABLED',
                                                         spatial_reference=sr)
        arcpy.AddField_management(connectors, 'FEATUREID', "LONG")
        arcpy.AddField_management(connectors, origin_feat_field, "LONG")
        arcpy.AddField_management(connectors, dist_orig_field, "FLOAT")
        arcpy.AddField_management(connectors, 'RPTDAYS', "FLOAT")
        arcpy.AddField_management(connectors, spatial_band_field, "FLOAT")
        arcpy.AddField_management(connectors, temporal_band_field, "FLOAT")

        # Insert connecting line features from the array of values
        fields = ['FEATUREID', origin_feat_field, dist_orig_field, 'RPTDAYS', spatial_band_field, temporal_band_field, 'SHAPE@']
        with arcpy.da.InsertCursor(connectors, fields) as rows:
            for new_line in new_lines:
                rows.insertRow(new_line)

        # Overall incident counts
        inc_cnt = 0
        orig_cnt = 0
        rpt_cnt = 0
        nrpt_cnt = 0

        # Classify & count incidents by type
        origins = list(set(oids) - set(nrids))  # origins excluding nr and r
        fields = ["OID@", incident_type_field, date_field, z_value_field]

        with arcpy.da.UpdateCursor(in_features, fields) as rows:
            for row in rows:
                inc_class = ''
                if row[0] in origins:
                    inc_class = 'O'
                    orig_cnt += 1
                elif row[0] in rids:
                    inc_class = 'R'
                    rpt_cnt += 1
                elif row[0] in nrids:
                    inc_class = 'NR'
                    nrpt_cnt += 1

                if inc_class:
                    row[1] = inc_class
                    inc_cnt += 1

                # calc z value if missing
                if not row[3]:
                    zval = row[2] - min_date
                    row[3] = zval.days

                rows.updateRow(row)

        # Get unit of feature class spatial reference system
        try:
            unit = units[sr.linearUnitName]
        except KeyError:
            unit = ''

        # Get half-life and half-distance
        test_distances = []
        half_distances = {}
        for sband in spatial_bands:
            test_distances.extend(all_distances[sband])
            test_distances.sort()
            half_distances[sband] = test_distances[int(len(test_distances)/2)]

        test_lives = []
        half_lives = {}
        for tband in temporal_bands:
            test_lives.extend(all_lives[tband])
            test_lives.sort()
            half_lives[tband] = test_lives[int(len(test_lives)/2)]

        # Build report content
        perc_o = "{:.1f}".format(100.0*float(orig_cnt)/inc_cnt)
        perc_nr = "{:.1f}".format(100.0*float(nrpt_cnt)/inc_cnt)
        perc_r = "{:.1f}".format(100.0*float(rpt_cnt)/inc_cnt)

        report_header = ('Repeat and Near Repeat Incident Summary\n'
                         'Created {}\n'.format(now_nice))

        data_info = ('Data Source: {}\n'
                     'Incident Date Range: {} - {}\n'.format(in_features, min_date, max_date))

        inc_type_report = ('Count and percentage of each type of incident\n'
                           ', Count, Percentage\n'
                           'All Incidents,{}, 100\n'
                           'Originators,{},{}\n'
                           'Near Repeats,{},{}\n'
                           'Repeats,{},{}\n'.format(inc_cnt,
                                                    orig_cnt, perc_o,
                                                    nrpt_cnt, perc_nr,
                                                    rpt_cnt, perc_r))
        console_type_rpt = ('Count and percentage of each type of incident\n'
                           '                  Count      Percentage\n'
                           'All Incidents   {:^10} {:^13}\n'
                           'Originators     {:^10} {:^13}\n'
                           'Near Repeats    {:^10} {:^13}\n'
                           'Repeats         {:^10} {:^13}\n'.format(inc_cnt, 100,
                                                                  orig_cnt, perc_o,
                                                                  nrpt_cnt, perc_nr,
                                                                  rpt_cnt, perc_r))

        half_lives_str = 'Estimated incident half-life\n'
        half_lives_str_console = 'Estimated incident half-life\n'
        for tband in temporal_bands:
            half_lives_str += '{} days temporal band, {:.1f} days\n'.format(tband, half_lives[tband])
            half_lives_str_console += '{} days temporal band: {:.1f} days\n'.format(tband, half_lives[tband])

        half_distance_str = 'Estimated incident half-distance\n'
        half_distance_str_console = 'Estimated incident half-distance\n'
        for sband in spatial_bands[1:]:
            half_distance_str += '{0} {1} spatial band, {2:.1f} {1}\n'.format(sband, unit, half_distances[sband])
            half_distance_str_console += '{0} {1} spatial band: {2:.1f} {1}\n'.format(sband, unit, half_distances[sband])

        temp_band_strs = ["<{} days".format(b) for b in temporal_bands]
        temporal_band_labels = ','.join(temp_band_strs)
        console_tband_labels = ' '.join(['{:^12}'.format(bnd) for bnd in temp_band_strs])

        counts_title = 'Number of Repeat and Near-Repeat incidents per spatial and temporal band\n'
        percent_title = 'Percentage of all incidents classified as Repeat or Near-Repeat and appearing in each spatial and temporal band\n'

        counts_header = ',{}\n'.format(temporal_band_labels)
        console_counts_header = '                 {}'.format(console_tband_labels)

        percent_header = ',{}\n'.format(temporal_band_labels)
        console_perc_header = '                 {}'.format(console_tband_labels)

        counts_table = ""
        percent_table = ""
        console_count = ""
        console_perc = ""

        row_sum = [0 for tband in temporal_bands]

        for sband in spatial_bands:

            # get temporal bands and their incident counts
            vals = band_counts[sband]

            # Get spatial band count in each temporal band
            # Sums include counts from smaller bands
            row_counts = [vals[tband] for tband in temporal_bands]
            try:
                row_sums = [sum(row_counts[0:i]) for i in xrange(1,len(row_counts)+1)]
            except:
                row_sums = [sum(row_counts[0:i]) for i in range(1,len(row_counts)+1)]

            row_sum = [x + y for (x, y) in zip(row_sums, row_sum)]
            row_perc = [100.0 * float(val)/inc_cnt for val in row_sum]

            # append counts & percentages to the table
            counts_table += '<{} {},{}\n'.format(sband, unit, ','.join([str(cnt) for cnt in row_sum]))
            console_count += '{:>16} {}\n'.format('<{} {}'.format(sband, unit), ' '.join(['{:^12}'.format(cnt) for cnt in row_sum]))
            percent_table += '<{} {},{}\n'.format(sband, unit, ','.join(["{:.1f}".format(prc) for prc in row_perc]))
            console_perc += '{:>16} {}\n'.format('<{} {}'.format(sband, unit), ' '.join(['{:^12}'.format("{:.1f}".format(prc)) for prc in row_perc]))

        # Write report
        reportname = path.join(report_location, "{}_{}.csv".format('Summary', now))
        with open(reportname, 'w') as report:

            report.write(report_header)
            report.write('\n')
            report.write(data_info)
            report.write('\n')
            report.write(half_distance_str)
            report.write('\n')
            report.write(half_lives_str)
            report.write('\n')
            report.write(inc_type_report)
            report.write('\n')
            report.write(counts_title)
            report.write(counts_header)
            report.write(counts_table)
            report.write('\n')
            report.write(percent_title)
            report.write(percent_header)
            report.write(percent_table)

        arcpy.SetParameterAsText(9, path.join(out_lines_dir, out_lines_name))
        arcpy.AddMessage("\nView incident summary report: {}\n".format(reportname))

        arcpy.AddMessage(report_header)
        arcpy.AddMessage('')
        arcpy.AddMessage(data_info)
        arcpy.AddMessage('')
        arcpy.AddMessage(half_distance_str_console)
        arcpy.AddMessage('')
        arcpy.AddMessage(half_lives_str_console)
        arcpy.AddMessage('')
        arcpy.AddMessage(console_type_rpt)
        arcpy.AddMessage('')
        arcpy.AddMessage(counts_title)
        arcpy.AddMessage(console_counts_header)
        arcpy.AddMessage(console_count)
        arcpy.AddMessage('')
        arcpy.AddMessage(percent_title)
        arcpy.AddMessage(console_perc_header)
        arcpy.AddMessage(console_perc)

    except arcpy.ExecuteError:
        # Get the tool error messages
        msgs = arcpy.GetMessages()
        arcpy.AddError(msgs)
        print(msgs)

    except:
        # Return  error messages for use in script tool or Python Window
        arcpy.AddError(str(sys.exc_info()[1]))

        # Print Python error messages for use in Python / Python Window
        print(str(sys.exc_info()[1]) + "\n")


if __name__ == '__main__':
    argv = tuple(arcpy.GetParameterAsText(i)
                 for i in range(arcpy.GetArgumentCount()))
    classify_incidents(*argv)
