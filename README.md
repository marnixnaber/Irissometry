# Irissometry

Thank you for downloading the irissometry toolbox!

The irissometry implementation does the following:
- Detects the eye's pupil in close-up videos using a starburst-like algorithm (see function detectPupil)
- Tracks points and calculate distances within the iris
- Spits out a matrix containing data about the pupil center coordinates, pupil size, iris feature distances, etc.
- Saves output matrix in .mat file and .csv file
- Plots some graphs with pupil radius and position data over time

For more info about input and output (io), enter "help irissometry" in the command.

See "example.m" for an example code. Just run it, select one or more videos, and observe the magic.

> Please cite our work in case you use our implementation!

![Example of irissometry output](https://github.com/marnixnaber/Irissometry/blob/main/images/irissometry.png)

## Reference:
Strauch, C., & Naber, M. (Submitted). Irissometry: effects of pupil size on iris elasticity measured with video-based feature tracking. Investigative Ophthalmology and Vision Sciences.

## Info for videos:
Make sure that the videos do not display black borders at the edge of the frame 
because this causes the pupil border detection to fail. This is how a video frame should look like:

![Example of a good video](https://github.com/marnixnaber/Irissometry/blob/main/images/goodVideoForIrissometry.png)

## Required software:
The code has been tested in MATLAB 2019b on a windows 10 machine. 
No guarantees can be provided for other MATLAB versions and operating platforms.

Please note that you need the computer vision toolbox and image processing toolbox to get the code to work.

Please contact Marnix Naber in case you cannot get the code to work.

## License:
This software can used for non-commercial purposes (see license file). Please contact marnixnaber@gmail.com for commercial licences.

## Acknowledgments:
Richard Brown (2007) for sharing the circle fit code.

Dongheng, L.,Winfield, D., and Parkhurst, D. J. Starburst:
A hybrid algorithm for video-based eye tracking combining
feature-based and model-based approaches, in 2005
IEEE Computer Society Conference on Computer Vision
and Pattern Recognition (CVPR'05) - Workshops, 2005,
pp. 79–79.

## Contact:
For questions, please contact marnixnaber@gmail.com
