# Testing Data
The bogus data includes randomly generated data with silly pictures for testing the viewer.

- All bogus databases contain occasional NaN or undefined values.
- Some databases will contain references to 0.png, which does not exist and should appropriately display the "broken image" icon where the image should be.
- Some databases contain references to abc.csv which is not a valid image and should display "could not display file: abc.csv" where the image should.
- bogus\_1 has an axis_order.csv file for testing axis orders.
- bogus\_1 also has 'query' and 'selection' parameters defined.
- big_bogus 2 through 4 are "very large" datasets (300+ rows). The scatterPlot and pcoord components should switch to their canvas versions when loading these.

## To Do:
- Add intentionally not-to-spec datasets to test error cases.