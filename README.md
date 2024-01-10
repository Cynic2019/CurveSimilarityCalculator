# Curve Similarity Calculator
A Python library for calculating the similarity of 1D, 2D, and 3D curves.

## Description
This project implements the calculation of similarity for 1D, 2D, and 3D curves, with results ranging from [0, 1].

- 1D refers to one-dimensional lines (arrays formed by a series of x coordinates) and utilizes the Levenshtein distance for similarity computation.

- 2D denotes two-dimensional curves (curves formed by a series of x, y coordinates). After Procrustes normalization and rotation of the curve, the Fréchet distance is used to compute similarity.

- 3D describes three-dimensional curves (curves formed by a series of x, y, z coordinates). After Procrustes normalization and matrix transformation of the curve, the Fréchet distance is used for similarity calculation.

## System Requirements
Python 3.5 or above and related dependency packages.

If your Python version is below 3.5, please modify the `curve_types.py` file by changing
```
Curve = List[Point]
```
to
```
Curve = list
```
and remove the corresponding dependency import: 
```
from typing import List
```
This is because `typing` is a standard library in Python 3.5 and later.

## Start the Project
1. Setup
   - Clone or download the project source code.
   - Import the project into your preferred IDE.
   - Install the required libraries.

2. Sample Data
   - For 1D curves, a simple template data is provided directly in the __main__ of curve_similarity_1d.py.
   - For 2D and 3D curves, sample data files are provided in the sample_data folder:
      - 2D curves use the 'Longitude' and 'Latitude' for comparison.
      - 3D curves use 'Longitude', 'Latitude', and 'Height' for comparison.
   - Upon setting up, you can run the project and expect the following similarity outputs:
     ```
     1D similarity: 0.6
     2D similarity: 0.8944347521625984
     3D similarity: 0.9140664484123516
     ```

3. Using Your Own Data

   For 1D curves, if your curve is lengthy, it's recommended to change the hard-coded str1 and str2 in the __main__ of curve_similarity_1d.py to read from a file.

   For 2D and 3D curves, organize your data into a format that pandas can read.  
   Then, in curve_similarity_2d.py and curve_similarity_3d.py, modify the file-reading method in __main__ and adjust the column names like 'Longitude', 'Latitude', and 'Height' accordingly.
   
## Steps and Plot

***3d***

1.Original Curve
![image](https://github.com/Cynic2019/CurveSimilarityCalculator/assets/37062530/aff7d9cb-8cad-4299-9ea3-bf1a2f1b2f8d)
2.Resample Curve
![image](https://github.com/Cynic2019/CurveSimilarityCalculator/assets/37062530/7fbe3d11-bcbf-481e-a9ba-eca93f5cea52)
3.Translate Curve
![image](https://github.com/Cynic2019/CurveSimilarityCalculator/assets/37062530/0c472d01-bd26-4c85-93d6-b412bd420335)
4.Zoom Curve
![image](https://github.com/Cynic2019/CurveSimilarityCalculator/assets/37062530/f33a01de-1a9f-40a2-8daf-6f17a72c972c)
5.Matrix-Transform Curve
![image](https://github.com/Cynic2019/CurveSimilarityCalculator/assets/37062530/d3acdc72-60f5-4a34-8070-a163ac9d983e)
The similarity between the two curves was calculated as: 0.9140664484123516

## Plot Code
Readers can use Matplotlib, Plotly, and Mayavi for plotting.

The diagram in my example was drawn using Plot, and the code is as follows:
```
def save_3D_curves_plotly(curve1, curve2, title, filename):
    curve1 = np.array(curve1)
    curve2 = np.array(curve2)

    # create a 3D curve chart
    fig = go.Figure()

    # add curve1
    fig.add_trace(go.Scatter3d(x=curve1[:, 0], y=curve1[:, 1], z=curve1[:, 2],
                                mode='lines',
                                line=dict(color='blue', width=5),
                                name='Curve 1'))

    # add curve2
    fig.add_trace(go.Scatter3d(x=curve2[:, 0], y=curve2[:, 1], z=curve2[:, 2],
                                mode='lines',
                                line=dict(color='red', width=5),
                                name='Curve 2'))

    print("pio.renderers.default: " + pio.renderers.default)

    # set layout properties
    fig.update_layout(title=title,
                      scene=dict(
                          xaxis_title='X',
                          yaxis_title='Y',
                          zaxis_title='Z'
                      ),
                      scene_camera = dict(
                          eye=dict(x=1.5, y=1.5, z=1.5)
                      )
                )

    # save as still image
    fig.write_image(filename + '.svg')

    # display graphics
    fig.show()
```
If you want to draw a graph of each step of curve change, you need to transfer the set of curve points after each step of curve change into the drawing function.

## License
This project is released under the MIT license.

## Contributing
Contributions are welcome! 

If you find any bugs or new requirements in the code, please feel free to contact me and develop together.

## Contact Information
If you have any questions about the project, feel free to contact me at:  
liyajing20@mails.tsinghua.edu.cn
