# 3D Object Tracking

- The code is for object tracking in 3D image volume sequences.

- The input should be a MAT file containing a 3D volume sequence. In each volume, background voxels are labeled as zero, and the voxels belonging to the same object are assign with the same but unique label.

- The output is similar, but the labels of each object in different frames will keep consistent.

## Reference
[Hanyi Yu, Fusheng Wang, Sung Bo Yoon, Robert Kauffman, Jens Wrammert, Adam Marcus, Jun Kong, “Non-Gaussian Models for Object Motion Analysis with Time-lapse Fluorescence Microscopy Images,” Modern Statistical Methods for Health Research, pp:15-41, Springer, 2021.](https://link.springer.com/chapter/10.1007/978-3-030-72437-5_2)

## License
This tool is available under the GNU General Public License (GPL) (https://www.gnu.org/licenses/gpl-3.0.en.html) and the LGPL (https://www.gnu.org/licenses/lgpl-3.0.en.html).
