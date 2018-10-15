This is a Python implementation of Dr. Sumner's Deformation Transfer for Triangle Meshes paper.
Note that the algorithm for finding correspondence is not implemented since we're using our own methods to do so.
Special thanks to Dr. Sumner himself for helping me out. 
Usage: Scroll down to the main function and change your input/output paths to your .obj files. 
The example included transfers a sad face to a laugh face. Note that in this case source and target base meshes are the same.
The program will only write positions of vertices to output. You need to paste vertex normal information and triangle indices from
the file vnandf (these are the same for any mesh included in this example though)
