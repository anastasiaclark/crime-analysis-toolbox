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
#              ArcGIS Advanced license required
#              Python 2.7
# author: ArcGIS Solutions
# contact: ArcGISTeamLocalGov@esri.com
# company: Esri
# ==================================================
# description: Classify incidents according to spatial and temporal proximity
#              to preceeding incidents. Generate a summary report of the
#              processed incidents
# ==================================================
# history:
# v1 05/05/2016
# ==================================================

import arcpy
from datetime import datetime as dt
from datetime import timedelta as td
from os import path

# Added field names
spatial_band_field = 'SPACEBAND'
temporal_band_field = 'TIMEBAND'
origin_feat_field = 'ORIGIN'
z_value_field = 'ZVALUE'
dist_orig_field = "DISTTOORIG"

units = {"Meter": {"short":"m",
                   "long": "Meters"},
         "Foot_US": {"short":"ft",
                     "long": "Feet"},
         "Foot": {"short":"ft",
                     "long": "Feet"}}


def reset_fields(fc):
    """Checks for required incident classification fields,
       and deletes/adds fields as necessary"""

    # Delete classification fields if they already exist in the dataset
    inc_fields = [f.name for f in arcpy.ListFields(fc)]

    delete_fields = []

    for field in[z_value_field]:
        if field in inc_fields:
            delete_fields.append(field)

    if delete_fields:
        arcpy.DeleteField_management(fc, delete_fields)

    # Add field for z value calculation
    arcpy.AddField_management(fc,
                              field_name=z_value_field,
                              field_type='LONG')


def calculate_band(value, bands):
    """Finds the first number in a list larger than a value"""
    for band in bands:
        if band > value:
            return band


def get_date_range(fc, field):
    """Returns the min and max values from a date field"""
    with arcpy.da.SearchCursor(fc, field) as rows:
        date_vals = [row[0] for row in rows if row[0]]

    date_vals = list(set(date_vals))
    date_vals.sort()

    # Range of incident dates
    return date_vals[0], date_vals[-1]


def find_feature_pair(fc, oid_field, date_field, ids, min_date):
    """Returns coordinates and attributes for vertices of line features"""

    links = []

    # Find the two features that are part of the connection
    where_clause = """{} in ({},{})""".format(oid_field, ids[0], ids[1])
    fields = [oid_field, date_field, z_value_field, 'SHAPE@X', 'SHAPE@Y']
    with arcpy.da.UpdateCursor(fc,
                               field_names=fields,
                               where_clause=where_clause) as cur_link:
        for feat in cur_link:
            # Calculate the z values of each incident in the pair
            if not feat[1]:
                return ['error', feat[0]]
            zval = feat[1] - min_date
            feat[2] = zval.days
            cur_link.updateRow(feat)
            links.append([feat[0], feat[1], feat[3], feat[4], feat[2]])

    return links


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

        # test for input features
        in_count = arcpy.GetCount_management(in_features)
        if int(in_count.getOutput(0)) <= 1:
            raise Exception("Insufficient features for processing")


        # Build sorted lists of unique band values
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

        # Get min and max date values
        min_date, max_date = get_date_range(in_features, date_field)

        # Get unit of feature class spatial reference system
        sr = arcpy.Describe(in_features).spatialReference
        try:
            unit_short = units[sr.linearUnitName]["short"]
            unit_long = units[sr.linearUnitName]["long"]
        except KeyError:
            unit_short = ''
            unit_long = ''

        # Keep track of origins and nrs
        oids = []
        nrids = []
        rids = []

        # Connecting line segments and table rows
        new_lines = []
        new_rows = []

        # Build empty dictionary to hold type tallies
        type_counts = {}
        for sband in spatial_bands:
            type_counts[sband] = {}
            for tband in temporal_bands:
                type_counts[sband][tband] = {'oids': [],
                                             'nrids': [],
                                             'rids': []}

        # Value lists for half life calculations
        all_distances = {}
        for sband in spatial_bands:
            all_distances[sband] = []

        all_lives = {}
        for tband in temporal_bands:
            all_lives[tband] = []

        found_connections = []

        # Build table of records within the max spatial band of another feature
        search = "{} {}".format(spatial_bands[-1], unit_long)
        near_table = arcpy.GenerateNearTable_analysis(in_features,
                                                      in_features,
                                                      search_radius=search,
                                                      closest='ALL',
                                                      method='GEODESIC')

        # Identify and process relevent near features
        fields = ['IN_FID', 'NEAR_FID', 'NEAR_DIST']
        ignored = []
        with arcpy.da.SearchCursor(near_table, field_names=fields) as nearrows:

            # Process each identified connection within the spatial bands
            for nearrow in nearrows:
                dist = nearrow[2]

                links = find_feature_pair(in_features,
                                          oidname,
                                          date_field,
                                          [nearrow[0], nearrow[1]],
                                          min_date)

                # Identify which feature is the oldest and id it as the source
                if links[0] == 'error':
                    if links[1] not in ignored:
                        arcpy.AddWarning('Ignored potential connection(s) with feature {0} because of invalid date'.format(links[1]))
                        print('Ignored potential connection(s) with feature {0} because of invalid date'.format(links[1]))
                        ignored.append(links[1])
                    continue

                elif links[0][1] > links[1][1]:
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

                # Identify the spatial & temporal bands
                link_found = False
                for sband in spatial_bands:
                    if dist <= sband:
                        for tband in temporal_bands:
                            if daydiff <= tband:
                                if not link_found:
                                    # track distances and lives for half measures
                                    all_distances[sband].append(dist)
                                    all_lives[tband].append(daydiff)
                                    incident_sband = sband
                                    incident_tband = tband

                                    link_found = True

                                # id classification
                                if oid not in type_counts[sband][tband]['oids']:
                                    type_counts[sband][tband]['oids'].append(oid)
                                if dist <= spatial_bands[0]:
                                    if fid not in type_counts[sband][tband]['rids']:
                                        type_counts[sband][tband]['rids'].append(fid)
                                elif fid not in type_counts[sband][tband]['nrids']:
                                    type_counts[sband][tband]['nrids'].append(fid)

                if link_found:
                    found_connections.append((oid, fid))

                    # create connecting line from x, y, z values of two pts
                    end = arcpy.Point(X=fx, Y=fy, Z=fz)
                    start = arcpy.Point(X=ox, Y=oy, Z=oz)
                    vertices = arcpy.Array([start, end])
                    feature = arcpy.Polyline(vertices, None, True, False)
                    new_lines.append([fid, oid, dist, daydiff, incident_sband,
                                      incident_tband, feature])

        # Delete near table
        arcpy.Delete_management(near_table)

        # Create feature class for connecting lines
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
        fields = ['FEATUREID', origin_feat_field, dist_orig_field, 'RPTDAYS',
                  spatial_band_field, temporal_band_field, 'SHAPE@']
        with arcpy.da.InsertCursor(connectors, fields) as rows:
            for new_line in new_lines:
                rows.insertRow(new_line)

        # Manage classification fields
        fieldnames = []
        for sband in spatial_bands:
            for tband in temporal_bands:
                fieldnames.append('s{}t{}'.format(int(sband), int(tband)))

        cur_fields = [f.name for f in arcpy.ListFields(in_features)]
        for fieldname in fieldnames:
            if fieldname in cur_fields:
                arcpy.DeleteField_management(in_features, fieldname)
            arcpy.AddField_management(in_features, fieldname,
                                      'TEXT', field_length=2)

        # Classify & count incidents by type
        for sband in spatial_bands:
            for tband in temporal_bands:
                band = type_counts[sband][tband]
                type_counts[sband][tband]['oids'] = [id for id in band['oids'] if id not in band['nrids'] and id not in band['rids']]
                type_counts[sband][tband]['nrids'] = [id for id in band['nrids'] if id not in band['rids']]

        fields = ["OID@", date_field, z_value_field]
        fields.extend(fieldnames)

        with arcpy.da.UpdateCursor(in_features, fields) as rows:
            inc_count = 0
            for row in rows:

                # calc z value if missing
                if not row[2]:
                    if not row[1]:
                        arcpy.AddWarning('Could not calculate z value for feature {} because of invalid date'.format(row[0]))
                        print('Could not calculate z value for feature {} because of invalid date'.format(row[0]))
                        continue

                    zval = row[1] - min_date
                    row[2] = zval.days
                inc_count += 1

                classifications = []

                for sband in spatial_bands:
                    for tband in temporal_bands:
                        if row[0] in type_counts[sband][tband]['nrids']:
                            classifications.append('NR')
                        elif row[0] in type_counts[sband][tband]['rids']:
                            classifications.append('R')
                        elif row[0] in type_counts[sband][tband]['oids']:
                            classifications.append('O')
                        else:
                            classifications.append('')
                row[3:] = classifications

                rows.updateRow(row)

        # Build empty dictionary to hold spatial and temporal band tallies
        band_counts = {}
        for sband in spatial_bands:
            band_counts[sband] = {}
            for tband in temporal_bands:
                band_counts[sband][tband] = 0

        for sband in spatial_bands:
            for tband in temporal_bands:
                if sband == spatial_bands[0]:
                    band_counts[sband][tband] = len(type_counts[sband][tband]['rids'])
                else:
                    band_counts[sband][tband] = len(type_counts[sband][tband]['nrids'])

        # Get half-life and half-distance
        test_distances = []
        half_distances = {}
        for sband in spatial_bands:
            test_distances.extend(all_distances[sband])
            test_distances.sort()
            if len(test_distances) > 0:
                half_distances[sband] = test_distances[int(len(test_distances)/2)]
            else:
                half_distances[sband] = 'Not Calculated'

        test_lives = []
        half_lives = {}
        for tband in temporal_bands:
            test_lives.extend(all_lives[tband])
            test_lives.sort()
            if len(test_lives) > 0:
                half_lives[tband] = test_lives[int(len(test_lives)/2)]
            else:
                half_lives[tband] = 'Not Calculated'

        # Build report content
        report_header = ('Repeat and Near Repeat Incident Summary\n'
                         'Created {}\n'.format(now_nice))

        data_info = ('Data Source: {}\n'
                     'Incident Date Range: {} - {}\n'
                     '# Incidents Processed: {}'.format(in_features,
                                                        min_date,
                                                        max_date,
                                                        inc_count))

        half_lives_str = '\nEstimated incident half-life\n'
        half_lives_str_console = 'Estimated incident half-life\n'
        for tband in temporal_bands:
            try:
                half_lives_str += '{} days temporal band, {:.1f} days\n'.format(tband, half_lives[tband])
                half_lives_str_console += '{} days temporal band: {:.1f} days\n'.format(tband, half_lives[tband])
            except ValueError:
                half_lives_str += '{} days temporal band, {} days\n'.format(tband, half_lives[tband])
                half_lives_str_console += '{} days temporal band: {} days\n'.format(tband, half_lives[tband])

        half_distance_str = 'Estimated incident half-distance\n'
        half_distance_str_console = 'Estimated incident half-distance\n'
        for sband in spatial_bands[1:]:
            try:
                half_distance_str += '{0} {1} spatial band, {2:.1f} {1}\n'.format(sband, unit_short, half_distances[sband])
                half_distance_str_console += '{0} {1} spatial band: {2:.1f} {1}\n'.format(sband, unit_short, half_distances[sband])
            except ValueError:
                half_distance_str += '{0} {1} spatial band, {2} {1}\n'.format(sband, unit_short, half_distances[sband])
                half_distance_str_console += '{0} {1} spatial band: {2} {1}\n'.format(sband, unit_short, half_distances[sband])

        temp_band_strs = ["<={} days".format(b) for b in temporal_bands]
        temporal_band_labels = ','.join(temp_band_strs)
        console_tband_labels = ' '.join(['{:^12}'.format(bnd) for bnd in temp_band_strs])

        counts_title = ('Number of Repeat and Near-Repeat incidents per '
                        'spatial and temporal band\n')
        percent_title = ('Percentage of all incidents classified as Repeat or '
                         'Near-Repeat and appearing in each spatial and '
                         'temporal band\n')

        counts_header = ',{}\n'.format(temporal_band_labels)
        console_counts_header = '{}{}'.format(' '*26, console_tband_labels)

        percent_header = ',{}\n'.format(temporal_band_labels)
        console_perc_header = '{}{}'.format(' '*26, console_tband_labels)

        counts_table = ""
        percent_table = ""
        console_count = ""
        console_perc = ""

        row_sum = [0 for tband in temporal_bands]

        for sband in spatial_bands:

            # get temporal bands and their incident counts
            vals = [band_counts[sband][tband] for tband in temporal_bands]

            # Get spatial band count in each temporal band
            # Sums include counts from smaller bands
            row_perc = [100.0 * float(val)/inc_count for val in vals]

            # append counts & percentages to the table
            if sband == spatial_bands[0]:
                counts_table += '<={} {},{}\n'.format(sband, unit_short, ','.join([str(cnt) for cnt in vals]))
                console_count += '{:>25} {}\n'.format('<={} {}'.format(sband, unit_short), ' '.join(['{:^12}'.format(cnt) for cnt in vals]))
                percent_table += '<={} {},{}\n'.format(sband, unit_short, ','.join(["{:.1f}".format(prc) for prc in row_perc]))
                console_perc += '{:>25} {}\n'.format('<={} {}'.format(sband, unit_short), ' '.join(['{:^12}'.format("{:.1f}".format(prc)) for prc in row_perc]))
            else:
                counts_table += '>{} to {} {},{}\n'.format(spatial_bands[0], sband, unit_short, ','.join([str(cnt) for cnt in vals]))
                console_count += '{:>25} {}\n'.format('>{} to {} {}'.format(spatial_bands[0], sband, unit_short), ' '.join(['{:^12}'.format(cnt) for cnt in vals]))
                percent_table += '>{} to {} {},{}\n'.format(spatial_bands[0], sband, unit_short, ','.join(["{:.1f}".format(prc) for prc in row_perc]))
                console_perc += '{:>25} {}\n'.format('>{} to {} {}'.format(spatial_bands[0], sband, unit_short), ' '.join(['{:^12}'.format("{:.1f}".format(prc)) for prc in row_perc]))

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
