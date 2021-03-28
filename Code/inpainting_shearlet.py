import os
import itertools
import sys
import concurrent.futures


import numpy as np
import imageio
import skimage.external.tifffile
import PIL
from scipy.ndimage import maximum_filter

import pyshearlab as psl


class ShearletInpainter:
    """
    Implementation of patch wise shearlet inpainting with cos² weighting
    Utilizes Iterative Thresholding according to
    https://www.math.tu-berlin.de/fileadmin/i26_fg-kutyniok/Kutyniok/Papers/InpaintAlgSPIE.pdf
    """
    def __init__(self, iterations, stop_factor, n_scales, patch_size=512, stride=384, max_workers=8, downsampling_factor=1):
        """
        Constructor
        :param iterations: number of iterations
        :type iterations: int
        :param stop_factor: stop factor for thresholding (0<stop_factor<1)
        :type stop_factor: float
        :param n_scales: n_scales for shearlet system, as in shearlab.org
        :type nscales: int
        :param patch_size: size of patch must be >= 512
        :type patch_size: int
        :param stride: stride, i.e., amount by which the window is shifted
        :type stride: int
        :param max_workers: maximum number of processes used for inpainting
        :type max_workers: int
        :param downsampling_factor: specify by which factor the image patches are downsampled before inpainting, must be a power of 2
        :type downsampling_factor: int
        """
        shearlet_system = psl.SLgetShearletSystem2D(1, int(patch_size/downsampling_factor), int(patch_size/downsampling_factor), n_scales)
        self.patch_size = patch_size
        self.stride = stride
        self.border_size = (self.patch_size - self.stride)
        self.max_workers = max_workers
        self.downsampling_factor = downsampling_factor
        self.stop_factor = stop_factor
        self.iterations = iterations
        self.shearlet_system = shearlet_system
  
    def inpaint_img(self, img, mask):
        """
        Inpaint a image by patchwise shearlet inpainting with cos² weighting
        :param img: img to be inpainted
        :type img: numpy array
        :param mask: mask that is 0 in void regions and 1 otherwise
        :type mask: numpy array
        :return: inpainted image
        """
        if img.shape != mask.shape:
            raise ValueError("img and mask must have same sizes")

        img, mask = self.pad_arrays(img, mask)
        inpainted_img = np.zeros(img.shape)                          
        coords = iter(self.calculate_coordinates())
                               
        with concurrent.futures.ProcessPoolExecutor(max_workers=self.max_workers) as executer:
            futures = []
            self.init_processes(executer, img, mask, coords, futures)
            
            while futures:
                # Wait for the next future to complete.
                done, futures = concurrent.futures.wait(
                futures, return_when=concurrent.futures.FIRST_COMPLETED
                )
                num_finished = len(done)
                self.finish_processes(inpainted_img, done)
   
                # Schedule the next set of futures.  We don't want more than max_workers futures
                # in the pool at a time, to keep memory consumption down.
                self.add_processes(executer, futures, img, mask, coords, num_finished)
            
        return inpainted_img[:self.array_shape_y, :self.array_shape_x]
    
    def calc_array_sizes(self, img):
        """
        calculate the sizes of the padded arrays
        """
        self.array_shape_y, self.array_shape_x = img.shape
        self.num_patch_x = int(self.array_shape_x/self.stride) + 1*(self.array_shape_x%self.stride!=0)
        self.num_patch_y = int(self.array_shape_y/self.stride) + 1*(self.array_shape_y%self.stride!=0)
        self.padded_array_shape_x = self.num_patch_x * self.stride + 2 * self.border_size
        self.padded_array_shape_y = self.num_patch_y * self.stride + 2 * self.border_size
    
    def pad_arrays(self, img, mask):
        """
        zero pad the arrays, that the sizes fit the patch-wise processing scheme
        """
        self.calc_array_sizes(img)
        img = zero_pad(img, self.padded_array_shape_x, self.padded_array_shape_y)
        mask = zero_pad(mask, self.padded_array_shape_x, self.padded_array_shape_y)
        return (img, mask)
    
    def calculate_coordinates(self):
        coords = [(y*self.stride, y*self.stride + self.patch_size, x*self.stride, x*self.stride + self.patch_size)
                  for x in range(self.num_patch_x)
                  for y in range(self.num_patch_y)
                 ]
        return coords
        
    def init_processes(self, executer, img, mask, coords, futures):
        for c in itertools.islice(coords, self.max_workers):
            img_patch, mask_patch = self.preprocess( img[c[0]:c[1], c[2]:c[3]], mask[c[0]:c[1], c[2]:c[3]] )
            future = executer.submit(self.inpaint_patch, img_patch, mask_patch, c)
            futures.append(future)
    
    def finish_processes(self, inpainted_img, done):
        for fut in done:
            inpainted_patch, coord = fut.result()
            inpainted_patch = self.postprocess(inpainted_patch)
            inpainted_patch = self.apply_cosine_weighting(inpainted_patch, coord)
            inpainted_img[coord[0]:coord[1], coord[2]:coord[3]] = inpainted_img[coord[0]:coord[1], coord[2]:coord[3]] + inpainted_patch
            del fut
    
    def add_processes(self, executer, futures, img, mask, coords, num_to_add):
        for c in itertools.islice(coords, num_to_add):
            img_patch, mask_patch = self.preprocess( img[c[0]:c[1], c[2]:c[3]], mask[c[0]:c[1], c[2]:c[3]] )
            future = executer.submit(self.inpaint_patch, img_patch, mask_patch, c)
            futures.add(future)

   
    def apply_cosine_weighting(self, inpainted_patch, coord):
        """
        apply the cos² weighting
        :param inpainted_patch: image patch after shearlet inpainting
        :type inpainted_patch: numpy array
        :param coord: coordinates of the patch
        :type coord: tuple of ints
        :return: image patch, weighted by cos² redundancy weights
        """
        if coord[0]!=0:
            inpainted_patch = inpainted_patch*self.calc_top_weighting()
        if coord[1]!=self.padded_array_shape_y-1:
            inpainted_patch = inpainted_patch*self.calc_bottom_weighting()
        if coord[2]!=0:
            inpainted_patch = inpainted_patch*self.calc_left_weighting()
        if coord[3]!=self.padded_array_shape_x-1:
            inpainted_patch = inpainted_patch*self.calc_right_weighting()
        return inpainted_patch
        
    def calc_left_weighting(self):
        """
        precompute the left border cos² weighting image
        :return: weighting window
        """
        left_weighting = np.ones((self.patch_size,)*2)
        sine = np.sin(np.arange(self.border_size)/self.border_size*np.pi/2)**2
        sines = np.vstack((sine,)*self.patch_size)
        left_weighting[:,:self.border_size] = sines
        return left_weighting
    
    def calc_right_weighting(self):
        """
        precompute the right border cos² weighting image
        :return: weighting window
        """
        right_weighting = np.ones((self.patch_size,)*2)
        cosine = np.cos(np.arange(self.border_size)/self.border_size*np.pi/2)**2
        cosines = np.vstack((cosine,)*self.patch_size)
        right_weighting[:,-self.border_size:] = cosines
        return right_weighting
        
    def calc_top_weighting(self):
        """
        precompute the top border cos² weighting image
        :return: weighting window
        """
        top_weighting = np.ones((self.patch_size,)*2)
        sine = np.sin(np.arange(self.border_size)/self.border_size*np.pi/2)**2
        sines = np.vstack((sine,)*self.patch_size).T
        top_weighting[:self.border_size,:] = sines
        return top_weighting
        
    def calc_bottom_weighting(self):
        """
        precompute the bottom border cos² weighting image
        :return: weighting window
        """
        bottom_weighting = np.ones((self.patch_size,)*2)
        cosine = np.cos(np.arange(self.border_size)/self.border_size*np.pi/2)**2
        cosines = np.vstack((cosine,)*self.patch_size).T
        bottom_weighting[-self.border_size:,:] = cosines
        return bottom_weighting
            
    def shearlet_inpainting(self, masked_img, mask):
        """
        Inpaint a patch using iterative thresholding and pyshearlab
        :param masked_img: image that has been masked : void regions set to 0
        :type masked_img: numpy array
        :param mask: mask that is 0 in void regions and 1 otherwise
        :type mask: numpy array
        :return: tuple of inpainted image patch and patch coordinates
        """
        normalized_coeffs = psl.SLnormalizeCoefficients2D(psl.SLsheardec2D(masked_img, self.shearlet_system), self.shearlet_system)
        delta = np.max(np.abs(normalized_coeffs))
        decay = self.stop_factor**(1/(max(1, self.iterations-1)))
        inpainted_img = np.zeros(masked_img.shape)    

        for i in range(self.iterations):
            res = mask * (masked_img-inpainted_img)
            coeffs = psl.SLsheardec2D(inpainted_img+res, self.shearlet_system)
            coeffs = coeffs*(np.abs(psl.SLnormalizeCoefficients2D(coeffs, self.shearlet_system))>delta)
            inpainted_img = psl.SLshearrec2D(coeffs,self.shearlet_system)
            delta=delta*decay
                                         
        return inpainted_img            
            
    def downsample_mask(self, mask):
        """
        downsample mask
        :param mask: mask that is 0 in void regions and 1 otherwise
        :type mask: numpy array
        :return: downsampled mask
        """
        target_size = int(mask.shape[0]/self.downsampling_factor)
        filtered_mask = maximum_filter(mask, size=(self.downsampling_factor , self.downsampling_factor ))
        x, y = mask.shape
        mask = PIL.Image.fromarray(filtered_mask)
        mask = mask.resize((target_size, target_size), PIL.Image.NEAREST)
        mask = np.asarray(mask)
        return mask
        
    def downsample_img(self, img_patch):
        """
        downsample image
        :param img_patch: image patch to be inpainted
        :type img: numpy array
        :return: downsampled image
        """
        target_size = int(img_patch.shape[0]/self.downsampling_factor)
        img_patch = PIL.Image.fromarray(img_patch)
        img_patch = img_patch.resize((target_size, target_size), PIL.Image.LANCZOS)
        img_patch = np.asarray(img_patch)
        return img_patch
    
    def upsample_img(self, inpainted_patch):
        """
        upsample image
        :param inpainted_patch: inpainted image patch to be upsampled
        :type img: numpy array
        :return: upsampled image
        """
        target_size = inpainted_patch.shape[0] * self.downsampling_factor
        inpainted_patch = PIL.Image.fromarray(inpainted_patch)
        inpainted_patch = img.resize((target_size, target_size), PIL.Image.LANCZOS)
        inpainted_patch = np.asarray(inpainted_patch)
        return inpainted_patch            
            
    def inpaint_patch(self, img_patch, mask_patch, coords=None):
        """
        Inpaint a patch using iterative thresholding and pyshearlab
        :param img_patch: image patch that has been masked : void regions set to 0
        :type img_patch: numpy array
        :param mask_patch: mask patch that is 0 in void regions and 1 otherwise
        :type mask_patch: numpy array
        :param coords: coordinates of the patch
        :type coords: tuple of ints
        :return: tuple of inpainted image patch and patch coordinates
        """
        # subtract median
        inpainted_patch = self.shearlet_inpainting(img_patch, mask_patch)
        return inpainted_patch, coords
                                       
    def preprocess(self, img_patch, mask_patch):
        """
        preprocess: downsample mask and masked_img, if specified
        :param img_patch: image patch to be preprocessed
        :type img_patch: numpy array
        :param mask_patch: mask patch to be preprocessed
        :type mask_patch: numpy array
        :return masked_img, mask
        """
        if self.downsampling_factor!=1:
            mask_patch = self.downsample_mask(mask_patch)
            img_patch = self.downsample_img(img_patch)
            return img_patch, mask_patch
        else:
            return img_patch, mask_patch
        
    def postprocess(self, inpainted_img):
        """
        postprocess image: upsample, if specified
        :param inpainted_img: inpainted image patch
        :type inpainted_img: numpy array
        :return: inpainted_img
        """
        if self.downsampling_factor!=1:
            inpainted_img = self.upsample_img(inpainted_img)
        return inpainted_img                                               

    
def zero_pad(img, output_size_x, output_size_y):
    """
    zero pad an image to given size
    :param img: image to be zero padded
    :param output_size_x: size in x dimension of zero padded image
    :param output_size_y: size in y dimension of zero padded image
    :return: zero padded image
    """
    img_zero_padded = np.zeros((output_size_y, output_size_x))
    img_zero_padded[:img.shape[0], :img.shape[1]] = img
    return img_zero_padded


def perform_patchwise_inpainting(input_file, inv_mask_file, ice_mask_file, output_file, iterations=400, stop_factor=0.001, n_scales=4, patch_size=512, stride=384, max_workers=4, downsampling_factor=1, apply_mask=True):
    """
    perform patchwise inpainting with cos² weighting
    :param input_file: path to input file
    :type input_file: str
    :param inv_mask_file: path to inv mask file
    :type inv_mask_file: str
    :param ice_mask_file: path to ice mask file
    :type ice_mask_file: str
    :param output_file: path to output file
    :type output_file: str
    :param iterations: number of iterations for iterative thresholding
    :type iterations: int
    :param stop_factor: stopfactor for iterative thresholding
    :type stop_factor: float
    :param n_scales: n_scales for shearlet construction,  for details shearlab.org
    :type n_scales: int
    :param patch_size: patch size for inpainting
    :type patch_size: int
    :param stride: stride (amount by which the window should be shifted for patch extraction)
    :type stride: int
    :param max_workers: maximum number of worker processes to be used for multiprocessing scheme
    :type max_workers: int
    :param downsampling_factor: specify by which factor the image patches are downsampled before inpainting, must be a power of 2
    :type downsampling_factor: int
    :return: None
    """
    print("Start Shearlet Inpainting (this may take a while) ...")
    print("Preprocessing ...")
    Inpainter = ShearletInpainter(iterations, stop_factor, n_scales, 
                                  patch_size, stride, max_workers, downsampling_factor)

    inv_mask = skimage.external.tifffile.imread(inv_mask_file)
    ice_mask = skimage.external.tifffile.imread(ice_mask_file)
    img = skimage.external.tifffile.imread(input_file)    

    # set the nan values to 0 and apply mask  (nan is encoded as -3.4e38)
    img = img * (img > -1e37).astype(int)
    void_mask = inv_mask < -1e37
    ice_mask = ice_mask > -1e37
    mask = np.logical_and(void_mask, ice_mask).astype(int)
    if apply_mask:
        img = img * (mask)
    
    print("Inpainting ...")
    result = Inpainter.inpaint_img(img, mask)
    
    print("Postprocessing...")
    result[mask == 1] = img[mask == 1]
    skimage.external.tifffile.imsave(output_file, result)
    skimage.external.tifffile.imsave("masked_input_debug.tif", img)

    print("Sucessfully Finished!")

    
    



