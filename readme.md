# Cinema Spec-D Parallel Coordinates Viewer

## Version 1.9
### Author: Cameron Tauxe

A Parallel Coordinates-based Viewer for Spec-D Cinema Databases

# Usage

## Installation
To install Cinema:Explorer, copy the **cinema.html** and **cinema/**
directory to the same directory as your data. In a simple example, the installed
files might look like this:

```
   cinema.html
   cinema/
   database_01.cdb
   database_02.cdb
   database_03.cdb
```

The **cinema/** directory is the general directory for all cinema applications,
each of which has an application and version directory. This is designed to
allow cinema applications to reside next to data, with a minimum of clutter from
Cinema. For example:

```
    cinema/explorer/1.9
        css/
        js/
        lib/
        license.md
        databases.json
```

The cinema databases can be anywhere below the directory that contains the
**cinema.html** and **cinema/** directory. So, you could also have this:

```
   cinema.html
   cinema/
   data/
       database_01.cdb
       database_02.cdb
       database_03.cdb
```

- To add your data to the viewer, edit **cinema/explorer/\<version\>/databases.json** to include references to your data.
    - You may override the default location and name of the **databases.json**
      by passing arguments in the URL. The file must be located at or below
      the level of the **cinema.html** file. For example:
```
    file:///Users/me/data/cinema.html?databases=somedirectory/somename.json
```

- **databases.json** is a list of objects where each object must contain *at least* a 'name' field and 'directory' field. 'name' is the name of the databases that will be shown in the viewer. 'directory' is the path to the '.cdb' directory for your database.
- Due to browser security features, most browsers will not read files if they are not in the same directory or in a subdirectory with **main.html**
- Optionally, each object in **databases.json** may contain a 'filter' field. This field is a regular expression to filter out dimensions in the database from being displayed on the chart in the viewer. Any dimension whose name matches with the regex will be filtered out.
- By default, any dimension beginning with FILE will be filtered, but this can be overridden with the 'filter' field.

## Viewing Data

- Start the viewer by opening **main.html** in a web browser.
- By default, the first database listed in **databases.json** will be loaded. To load a different database, select it in the "Database Select" drop-down at the top of the page and press "Load." (For large databases (300+ rows) a warning will be displayed).
- The viewer is split into two sections, the top section contains the Parallel Coordinates chart (called the pcoord area for brevity) and the bottom section contains a component for viewing the data, either a spread of images, or a scatter plot depending on the selected tab (this is called the view area). The two sections can be resized by click-and-dragging on the black bar between them.

## Using the Paralell Coordinates Chart

- Each vertical axis on the chart represents one dimension in the database. Each horizontal (gray) line represents a row/data point. It crosses through each axis at the row's value for that dimension.
- If a data point contains a NaN value, then the line will pass through the red area below the axis marked 'NaN.' If a data point has an undefined value, then the line will not be drawn through that axis.
- Moving the mouse over a line will display a readout of that data point in the corner of the screen.
- The axes can be re-arranged by click-and-dragging on their titles. If the database you loaded contains a valid axis_order.csv file, the top of the page will include drop-downs to select from preset orderings of the axes.
- You can filter data by click-and-dragging vertically along an axis. This will create a selection and the chart will only show data that passes through the selected region on the axis. You can filter the selection further by creating selections on other axes. Remove a selection by clicking on the axis outside of the selected area.
- The component in the view area will automatically update to show only the selected data on the chart.

## Querying the Parallel Coordinates Chart

- Click the '>' button on the left-side of the page by the chart to open the query panel. Click the '<' button to close it again.
- The query panel contains a slider for every numeric (Integer or Float) dimension in the database. They can be used to select an arbitrary value along that dimension.
- The checkbox next to each slider indicates whether or not to include that dimension in your query.
- When you define an arbitrary data point using the sliders and checkboxes, a red dotted line will be drawn on the chart to represent the data you defined. Two light pink lines on either side of the red line give a rough approximation of the range of results that your query will return.
- The threshold for the query (how close a data point has to be to one you defined to be included) can be changed with the "Threshold" input at the top of the panel. Notice that the two pink lines on the chart move closer to or further away from the red line as you change the threshold.
- A query can be performed by pressing the "Find Similar" button at the top of the panel. If any data was found within the threshold, the selections on the chart will be automatically changed to encapsulate all the data found.

## Using the Image Spread View

- With the "Image Spread" tab selected, the view area contains the images for all the currently selected data.
- Each box represents a selected data point and contains an image for every 'FILE' dimension in the database. If the dimension's filetype is a valid image (GIF,PNG or JPEG), then the image will be displayed. Otherwise, a warning will be shown saying that the file could not be displayed.
- Moving the mouse over a box will display a readout of that data point in the corner of the screen as well as highlight the corresponding data point in the Parallel Coordinates Chart.
- You can click on an image to view the full-size version. With the full-size version open, click anywhere on the screen to close it again.
- The thumbnail size of the images can be changed with the "Image Size" slider in the top-right corner of the panel.
- You can re-arrange the order of the data displayed in the view by using the "Sort By" drop-down to select a dimension to sort the data by. By default, data is sorted in ascending order, but it can be changed to descending by checking the "Reverse Sort Order" checkbox.
- If the data shown extends beyond the screen, you can scroll to view the rest of the page of data.
- If there is too much data to display on a single page, a window will appear at the bottom of the panel to allow you to browse through the different pages. The number of data points shown on each page can be changed with the "Results Per Page" drop-down menu.

## Using the Scatter Plot View

- With the "Scatter Plot" tab selected, the view area contains a 2D Scatter Plot with points for all the currently selected data.
- Moving the mouse over a dot will display a readout of that data point in the corner of the screen as well as highlight the corresponding data point in the Parallel Coordinates Chart (and vice versa).
- The two dimensions used as the x and y axes in the chart can be changed with the drop-down menus on the left-side and bottom of the panel.
- If data contains NaN or undefined values in at least one of the two dimensions selected, the data will not be plotted and a warning will be displayed in the lower-right corner of the panel.

# Testing

If you make changes to either the viewer or the underlying Cinema Components library, please follow these steps and ensure that everything still works as expected.

- Load a database.
- Load a database in Safari (It's "special")
- Load a database that contains NaN and undefined values.
- Load a database that is not to-spec and confirm that an error appears.
- Load a database with an axis_order.csv file.
- Load a database with an axis_order.csv file that is not to-spec and confirm that the option to select axis ordering does not appear as well as a warning in the developer console.
- Load a database while the Scatter Plot tab is selected (as opposed to the default Image Spread)
- Load a large database (300+ rows) and confirm that the pcoord and scatterPlot components switch to their respective canvas-based versions.
- Load a database while the query panel is open
- Load a database while the smoothLines checkbox is unchecked.
- Make a query with both the SVG and Canvas versions of the pcoord component.
- Select and highlight data with both the SVG and Canvas versions of the pcoord and scatterPlot components.
- Toggle "Smooth Lines" on and off with both the SVG and Canvas versions of the pcoord component.
- Change the selection while not on the first page of results in the image spread.
- Mess with all the settings in the image spread component and ensure that sorting and pagination are accurate.

# Changelog
## V1.9
- Update to v2.6.1 of Cinema Components Library
- Can specify 'query' and 'selection' parameters in Databases.json
- Alternative json to Databases.json can be specified with HTTP paramters
## V1.8
- Now using v2.4 of Cinema Components library
- Scatter Plot will remember which two dimensions you were viewing when changing tabs
## V1.7
- Now using v2.3 of Cinema Components library
- Added Scatter Plot viewing mode
## V1.6
- Now using v2.2 of Cinema Components library
- Added support for axis_order.csv files in databases
## V1.5
- Significant refactoring. Now used v2.0 of Cinema Components library
- ImageSpread and Query Panel have be re-implemented as components in library
## V1.4.3
- Page selection widget changes to a more compact form when there are many pages of results
- Added "Reverse Sort Order" option
- Fixed results not being sorted correctly when changing the database being viewed
## V1.4.2
- Update to use v1.4.02 of Parallel Coordinates Component
## V1.4.1
- Now displays all FILE fields for each result
- Fixed sorting not working properly when undefined or NaN values were involved
- Update to use v1.4.01 of Parallel Coordinates Component
## V1.4
- Compatible with v1.1 of Cinema SpecD and v1.4 of Parallel Coordinates Component
## V1.3.1
- Added toggle for smoothing lines on chart
## V1.3
- Up to Spec with Cinema Spec-D V1.0. (Can read both Type I and Type II)
## V1.2
- Updated to use version 1.3.1 of Parallel Coordinates Component
- Added feature to allow user to define a custom result, and query for results similiar to it.
- Filetype and dimension filter are defined in databases.json
## V1.1
- Updated to use version 1.2 of Parallel Cooridinates Component
## V1.0
- Initial
