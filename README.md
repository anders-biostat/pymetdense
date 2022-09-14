# pymetdense
Python package that implements a python parser for [.metdense files](https://github.com/anders-biostat/metdense).

## Usage
Install with pip
```bash
$pip install pymetdense
```
and import the relevant class:
```python
from pymetdense import MetDenseFile
```

A metdense file can then be instantiated with:
```python
file = MetDenseFile("path/to/file.metdense")
```

Some helpers:
```python
first_chromosome_number_of_positions = file.get_chromosome_length("1")
cell_names = file.get_cell_names()
```

If the file contains more cells than you want to load, you can set a mask that is applied to all data output. For example, if you have a list of desired cells `cells`, you can set
```python
cell_names = file.get_cell_names()
mask = [cell_names.index(cell) for cell in cells]
file.set_cell_mask(mask)
```
and the output will then be the data of `cells` in the correct order.

The positions and data within a region (inclusive for upper and lower bound) on a chromosome can be loaded using
```python
chromosome = "1"
lower, upper = 1e7, 1.1e7
positions = file.read_positions(chromosome, lower, upper)
data = file.read_data(chromosome, lower, upper)
```
where `positions` is an array of all the available positions within [lower, upper] on the chromosome, and `data` is an array of shape (positions, cells) that contains the data from the file.

You can also load the positions and data for multiple ranges at once, using
```python
chromosome = "1"
ranges = [[1e7, 1.1e7],[1.2e7,1.3e7]]
data, pos = file.read_data_positions(chromosome, ranges, mode='block')
```
where `data` and `pos` now are lists of the same arrays as above for each region. Setting `mode='block'` loads all data within the lowest and highest position across all ranges at once and filteres the relevant positions before further processing, which will typically be faster but has a higher memory requirement. `mode='separate'` simply runs a loop over the single range functions for all ranges.
