This code is used for data analysis for the paper  
"Functional hyperemia drives fluid exchange in the paravascular space" by Kedarasetti et. al. 
https://doi.org/10.1186/s12987-020-00214-3

All the code in this repository works with Matlab® 2019 (Mathworks Inc.)

To analyze the image files and calculate the displacements run the file "Thy1_analysis_summary.m"
The information about the datasets should be stored in an excelsheet with the following columns:
Animal ID; Day ID ;File Number;Pixel Scale

All the data including the excel files is available upon request. Please email raviteja(_at_)psu(dot)edu (or) pjd17(_at_)psu(dot)edu

The resulting displacements calculated from running this file are manually verified using ImageJ and the locations are added to the excel file. 

The displacement_summary.m file is run to calculate the mean and standard deviation of normalized impulse responses (Fig 5n)

The files and copyright information is as follows:

"dftregistration"			Copyright (c) 2007, Manuel Guizar 
"displacement_summary"		Copyright (c) 2019 Ravi Kedarasetti
"inset.m" 				Copyright (c) 2010, Moshe Lindner 
"IR_analytic.m"			Copyright (c) 2013 Bingxing Huo
"read_tiff_downsample.m"		Copyright (c) 2019 Ravi Kedarasetti
"register_filter_unleak.m"		Copyright (c) 2019 Ravi Kedarasetti
"Thy1_analysis_summary.m"		Copyright (c) 2019 Ravi Kedarasetti
"thy1_displacement_fields_iterative.m"	Copyright (c) 2019 Ravi Kedarasetti
"velocity_proc2.m"			Copyright (c) 2013 Bingxing Huo
"vessel_diameter_line.m"		Copyright (c) 2019 Ravi Kedarasetti
"wave_filt.m"			Copyright (c) 2019 Ravi Kedarasetti

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
