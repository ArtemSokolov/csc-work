import tifffile
import zarr
import napari
import yaml

import dask.array as da
import pandas as pd

def single_channel_pyramid(tiff_path, channel):
    print('Importing channel ' + str(channel) + ' from ' + str(tiff_path))
    tiff = tifffile.TiffFile(tiff_path)

    pyramid = [
        zarr.open(s[channel].aszarr())
        for s in tiff.series[0].levels
        ]

    pyramid = [
        da.from_zarr(z)
        for z in pyramid
        ]

    return pyramid

# import markers.csv
markers = pd.read_csv('data/markers.csv')

# import image contrast settings
with open('data/contrast_limits.yml') as f:
    contrast = yaml.safe_load(f)

# import cell coordinates
print('Loading cell coordinates')
xy = pd.read_csv('data/WD-76845-097-ij_subtracted_50_qc.csv')
xy = xy[['CellID','X_centroid','Y_centroid']]

# import cluster assignments
cluster = pd.read_csv('data/hdbscan.csv').merge(xy, on='CellID')

# add DNA1 channel to Napari image viewer
file_path = 'data/WD-76845-097.ome.tif'
dna = single_channel_pyramid(file_path, 0)
viewer = napari.view_image(
     dna, rgb=False, blending='additive',
     colormap='gray', visible=True, opacity=0.2,
     name='DNA1', contrast_limits=[152.0, 49353.0]
     )

# add segmentation outlines
seg_path = 'data/cellRingOutlines.tif'
seg = single_channel_pyramid(seg_path, 0)
viewer.add_image(
    seg, rgb=False, blending='additive',
    colormap='red', visible=True, opacity=0.3,
    name='segmentation'
    )

abx_colors = {'anti_CD3': 'bop blue', 'anti_CD45RO': 'bop orange',
              'CD4_488': 'bop purple', 'CD45_PE': 'bop orange',
              'PD1_647': 'bop blue', 'CD20_488': 'bop orange',
              'CD68_555': 'bop blue', 'CD8a_660': 'bop orange',
              'CD163_488': 'bop blue', 'FOXP3_570': 'bop orange',
              'Keratin_570': 'bop blue', 'aSMA_660': 'bop orange',
              'PCNA_488': 'bop orange', 'Desmin_555': 'bop blue',
              'Ecad_488': 'bop purple', 'PDL1_647': 'bop orange',
              'Vimentin_555': 'bop blue'}

abx_channels = ['Ecad_488_cellRingMask', 'CD8a_660_cytoRingMask',
                'anti_CD3_cytoRingMask', 'Keratin_570_cellRingMask',
                'PCNA_488_nucleiRingMask']

# loop over antibodies of interest and add them to Napari image viewer
for ch in abx_channels:
    ch = ch.rsplit('_', 1)[0]
    channel_number = markers['channel_number'][markers['marker_name'] == ch]
    img = single_channel_pyramid(
        file_path,
        channel=(channel_number.item() - 1)
        )
    viewer.add_image(
        img, rgb=False, blending='additive',
        colormap=abx_colors[ch], visible=False,
        name=ch, contrast_limits=(contrast[ch][0], contrast[ch][1])
        )

# highlight cluster of interest
for cl in [5,10,17]:
    centroids = cluster[cluster.hdbscan == cl][['Y_centroid', 'X_centroid']]
    viewer.add_points(
            centroids, name=f'cluster{cl}',
            face_color='white',
            edge_color='white',
            edge_width=0.0, size=4.0,
            visible=False
            )

napari.run()
