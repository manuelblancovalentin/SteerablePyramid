# SteerablePyramid
Steerable Pyramid Builder, Visualizer &amp; Texture Synthesizer for MATLAB

I created this code for my MSc. Thesis in Oil Reservoir Image Data 
 processing as I needed to understand how the workflow proposed by 
 Portilla et al worked. Even though Portilla's original code (see here: 
 http://www.cns.nyu.edu/~eero/steerpyr/) works great, it was very 
 difficult to fully comprehend each step in the different processes of 
 building the steerable pyramid, characterizing the textures and 
 synthesizing them, as almost no comments nor intelligible variable names 
 were used on the code. 
 I studied the code for several weeks and implemented all processes on my 
 own. I commented all steps and tried to include as many references to 
 the original paper as possible; so If anyone wants/needs to 
 understand how this code works, it'll be a lil' bit easier. 
 I have implemented ALL functions from scratch and I only use 1 of 
 Portilla's original code (expand function). I have also included another 
 function ('dispPyramid') which displays the steerable pyramid in a more 
 friendly way (giving the impression of a real pyramid), so that if you 
 need to display the pyramid in any paper or work, you can simply use 
 this function almost out the box. 
 I have tested my code and Portilla's original code using the same input 
 images and same initial conditions (starting white noise) and they 
 provide the EXACT same results (I checked this value per value). 
 This code can be used mainly for three purposes: 
  
  1) Build a steerable pyramid, given an input image, and a number of 
  desired scales and orientations. A steerable pyramid is a technique 
  that uses a recursive filtering workflow to obtain the different 
  attributes (texturally speaking) that an image may have at different 
  scales and orientations. This workflow basically consists on taking the 
  input image, filtering it with a highpass filter, a series of 
  different-oriented 2D band-pass filters and a low-pass filter. The 
  low-pass component (what's left-over) is then downscaled using a factor 
  of 2 and all the process is repeated again (High-pass + band-pass + 
  LP). By doing so we are extracting textural information about 
  orientation of the features in our image at different scales. 
  
  2) Extract Features that characterize the texture on the input image at 
  different scales and orientations (based on the steerable pyramid built 
  previously). These feature are those described on Portilla's paper too. 
  As shown on that paper they seem to be sufficient to characterize 
  several types of textures well (indeed they seem to do its job so well 
  that we can actually synthesize almost identical textures by using 
  them, as I will explain on the following entry). On the other side, 
  Portilla and Simoncelli asure they are universal parameters (which 
  roughly means that you can use it in a bunch of different applications 
  without having to implement much changes). 
  
  3) Texture Synthesis. As described in Portilla's paper, we can use the 
  textural features extracted from the steerable pyramid (used to 
  characterize it) to create a totally artificial and synthetic image 
  that will have the same textural characteristics as the original one. 
  On the other side, it is also possible to extrapolate texture, creating a 
  Mask that will preserve some part of the original texture, and synthesize the rest of it.
