# Software de Observaciones Sintéticas SOS

Developed by SPA group at INAOE for molecular cloud analysis
https://inaoep.mx/~astropol/

SOS is a package of tools to easily display and analize the physical propierties of a synthesized molecular cloud. 
It includes:
- Mass calculation throguh virial, LTE and X-factor methods
- Estimation of column density H2, 13CO
- Map binning functions
- Line fitting: automatic and manual modes 
- Line integrated profiles of a map section
- Display moment zero maps
- Display binning maps of: mass, column density and line profiles
- Display molecular lines
- Data base of the observed molecular clouds parameters (so far only TaurusMC)

## Requirements 📋

- Python 3.6 or later
- Astropy >= 4.0
- Scipy >= 1.5.1

## Install

Clone or download-unzip this GitHub repository. Change the directory to the sos path and import it

```python
cd sos_path/
import sos
```
## Usage example

An example to fit lines and calculate physical parameters of a given molecular cloud.

### Full molecular cloud

First, specify the paths for the molecular cube data of the molecular cloud (eg 'TaurusMC'). 

```python
import sos

db = sos.mc_db_tools.mc_db(sos.DB_PATH)
db.get_mc_params('TaurusMC')

db.path13CO = '/home/marcial/Documents/SPA/sos/data_cube/filament_13co_YZ_Jansky-per-px.fits'
db.path12CO = '/home/marcial/Documents/SPA/sos/data_cube/filament_12co_YZ_Jansky-per-px.fits'
db.pathC18O = '/home/marcial/Documents/SPA/sos/data_cube/filament_c18o_YZ_Jansky-per-px.fits'
db.save_mc_db()
```

The above just update the paths of the molecular cloud object in the database. It just needs to be run once, because the path is already updated unless the path of any molecule has changed.

Now, create the molecular cloud object (mc), using their ID name (eg 'TaurusMC'):

```python
mc = sos.mc('TaurusMC')
```
To display all the ID molecular clouds of the database:

```python
db_objs = sos.mc_db(sos.DB_PATH) 
db_objs.mc_db.keys()
```

Then, create the integrated velocity of the whole map for all the molecular lines available

```python
mc.get_map_vels()
```

Fit the line for the molecules

```python
mc.line_fit('13CO')
mc.line_fit('12CO')
```
Both molecules are needed in order to get most of the physical propierties

Finally, calculate the physical parameters: masses (three different methods) and column densities:

```python
mc.get_gral_params()
```

Display the results:

```python
mc.summary()
```

### Binned molecular cloud

Divide the molecule map into nbins 

```python
mc.binning_mol('13CO', nbins=16)
mc.binning_mol('12CO', 16)
```

Fit the line for each bin

```python
mc.line_fit_binning()
```

Calculate the physical parameters for each bin as:

```python
mc.get_bins_params()
```

#### Plot the binned resuls

First create the plotter object

```python
# Get the binned data results
mc_bins = mc.mc_binning

# Get the M0 map of one molecule
m0_data, m0_header = mc.M0('13CO', save=False)

# Create plotter object
plt = sos.mc_plotter(mc_bins, m0_data, m0_header)
```

Plot the lines over the M0 map, defining the molecule line

```python
plt.plot_m0_line('13CO')
```

### Or...

Run the example.py file, which contains the code lines above:

```python
import sos

run example.py
```

## Wiki 📖

Ongoing...

## License 📄

This project is licensed under the mit license. [LICENSE.md](LICENSE.md) for more details.

