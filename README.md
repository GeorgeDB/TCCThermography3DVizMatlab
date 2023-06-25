<h1>Multimodal Fusion for 3D Visualization of Defects/Cracks Using Infrared Themography</h1>
<h2>Abstract</h2>
<p>Along with the many advances in various fields of engineering, the demand for quality control of parts has increased. This need is directly related to reducing the risks of critical failures, unexpected stops in production lines and, consequently, greater safety, depending 
on the case analyzed. In order to make this control possible, some analysis techniques havebeen developed, from direct tests on the part to the quality of the final product.</p>
<p>Infrared thermography (IRT) has stood out from other failure analysis techniques because it can provide internal (subsurface) structural information about the part while being nondestructive and fast. However, its result is a 2D image and can be difficult to interpret. To
facilitate the interpretation of the results, the present work aimed to present a multimodal fusion analysis method, based on the thermographic image and on a 3D model of the part.</p>
<p>Therefore, the 3D modeling of the part in a software, the projection of its points on the plane, according to the referential of the camera used in the thermographic test, and
subsequent correlation between the projection and the thermography data were performed. In this way, it was concluded that multimodal fusion is an interesting technique to be
applied and that it can be of great help for the structural analysis of the part.</p>

<h2>Idea</h2>
<p>Having the infrared images of the objects, the details of how the test was performed and a CAD model of the object, map all the pixels to one point is space, in a way that with the pixel
color, the time frame and the details of the test, it will be possible to determine the depth of the specific 2D pixel in the 3D world. In that way, with this depth and the color-temperature
direct relation, it is possible to have a 3D model with each point having a distinct temperature. Therefore, it will be possible to determine if there is any crack or subsurface defect based 
on the temperature (if there is none, the heat flow should be "normal", depending on the material).</p>

<h2>Data</h2>
<ul>
  <li>Peça.STL (object) → object CAD model saved as STL (stereolithography) file, commonly used for 3D printing and CAD exporting/importing</li>
  <li>finalCode.m → final MATLAB code used to import the model, the images and rebuild the 3D visualization after mapping the infrared image pixels into the model points</li>
  <li>modernPosit.m → code used to determine points to be used as reference for relation between object and image</li>
  <li>refinePDEMmesh.m → code used to create mesh from STL model</li>
</ul>

<h2>Results Examples</h2>
<h3>Reference Points Used to Make Relation Between Model and Image</h3>
<img src="https://github.com/GeorgeDB/TCCThermography3DVizMatlab/assets/124641422/e4ada1db-e123-4c63-9df5-0f8a7ab40db2.png">
<h3>Mesh Nodes Printed Over Image</h3>
<p>Each black point correspond to a mesh node (model point) printed over the 2D image</p>
<img src="https://github.com/GeorgeDB/TCCThermography3DVizMatlab/assets/124641422/8578a7e3-686a-4717-9b41-38b7b55832e1.png">
<h3>Final Rebuild Model With Temperature Information</h3>
<img src="https://github.com/GeorgeDB/TCCThermography3DVizMatlab/assets/124641422/c64fcfc1-c55d-4909-bc43-38da6c2bbbb3.png">
