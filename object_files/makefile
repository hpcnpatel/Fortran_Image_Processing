#make file by Nisarg

dir= ../source_files/
obj=o

#----------------------------------------------------------------
#----------------------------------------------------------------

COM=ifort
PARALLEL=-qopenmp
#OPT=-O3
#FLAGS= -Wall -fbacktrace -fbounds-check -pg -g
FLAGS= -xHOST -O3 -g -traceback -check bounds

#FLAGS is a varible which can be called from cmd line as: "makefile FLAGS=-O3"
#f90=gfortran -fopenmp -Wall -fbacktrace -fbounds-check -std=f2003

COMPILER=$(COM) $(PARALLEL) $(FLAGS)

#----------------------------------------------------------------
#----------------------------------------------------------------

#all_objects =factorial.$(obj) \
#mod_vtkio.$(obj)\
#LIB_VTK_IO.$(obj)\

all_objects =factorial.$(obj) \
constants_module.$(obj)\
fftpack_module.$(obj)\
fast_fourier_transform_module.$(obj)\
general_fun_sub.$(obj)\
variable_module.$(obj)\
utilities_module.$(obj)\
ghost.$(obj)\
sort_module.$(obj)\
mod_gnuplot.$(obj)\
mod_vandermonde.$(obj)\
read_a_file.$(obj)\
binary_reader.$(obj)\
vtk_reader.$(obj)\
vtk_writter.$(obj)\
raw_reader.$(obj)\
mod_interpolate_ctdata.$(obj)\
filter.$(obj)\
QR_module_short.$(obj)\
QR_module_short_QR_patch.$(obj)\
mod_grad_of_image.$(obj)\
mod_median.$(obj)\
mod_qr.$(obj)\
mod_gauss_filter.$(obj)\
mod_surface.$(obj)\
image_process.$(obj)

#----------------------------------------------------------------
#----------------------------------------------------------------

inter.exe : $(all_objects)
	$(COMPILER) $(all_objects) -o ../inter

#----------------------------------------------------------------
#----------------------------------------------------------------

factorial.$(obj) :$(dir)factorial.f90
	$(COMPILER) -c  $(?)

variable_declaration.$(obj) :$(dir)variable_declaration.f90
	$(COMPILER) -c  $(dir)variable_declaration.f90

general_fun_sub.$(obj) :$(dir)general_fun_sub.f90
	$(COMPILER) -c  $(dir)general_fun_sub.f90

LIB_VTK_IO.$(obj) :$(dir)LIB_VTK_IO.f90
	$(COMPILER) -c  $(dir)LIB_VTK_IO.f90

ghost.$(obj) :$(dir)ghost.f90
	$(COMPILER) -c  $(dir)ghost.f90

mod_surface.$(obj):$(dir)mod_surface.f90
	$(COMPILER) -c $(dir)mod_surface.f90

vtk_writter.$(obj) : $(dir)vtk_writter.f90
	$(COMPILER) -c $(dir)vtk_writter.f90

binary_reader.$(obj) : $(dir)binary_reader.f90
	$(COMPILER) -c $(dir)binary_reader.f90

mod_type.$(obj) : $(dir)mod_type.f90
	$(COMPILER) -c $(dir)mod_type.f90

constants_module.$(obj) : $(dir)constants_module.f90
	$(COMPILER) -c $(dir)constants_module.f90

variable_module.$(obj) : $(dir)variable_module.f90
	$(COMPILER) -c $(dir)variable_module.f90

sort_module.$(obj) : $(dir)sort_module.f90
	$(COMPILER) -c $(dir)sort_module.f90

utilities_module.$(obj) :$(dir)utilities_module.f90
	$(COMPILER) -c $(dir)utilities_module.f90

mod_gnuplot.$(obj) :$(dir)mod_gnuplot.f90
	$(COMPILER) -c $(dir)mod_gnuplot.f90

mod_grad_of_image.$(obj) :$(dir)mod_grad_of_image.f90
	$(COMPILER) -c $(dir)mod_grad_of_image.f90

QR_module_short.$(obj) : $(dir)QR_module_short.f90
	$(COMPILER) -c $(dir)QR_module_short.f90

QR_module_short_QR_patch.$(obj) : $(dir)QR_module_short_QR_patch.f90
	$(COMPILER) -c $(dir)QR_module_short_QR_patch.f90

#mod_read_ct_data.mod : mod_read_ct_data.$(obj) $(dir)mod_read_ct_data.f90
#mod_read_ct_data.$(obj) : $(dir)mod_read_ct_data.f90
#	$(COMPILER) -c $(dir)mod_read_ct_data.f90

raw_reader.$(obj) : $(dir)raw_reader.f90
	$(COMPILER) -c $(dir)raw_reader.f90

read_a_file.$(obj) : $(dir)read_a_file.f90
	$(COMPILER) -c $(dir)read_a_file.f90

vtk_reader.$(obj) : $(dir)vtk_reader.f90
	$(COMPILER) -c $(dir)vtk_reader.f90

mod_interpolate_ctdata.$(obj) :   $(dir)mod_interpolate_ctdata.f90
	$(COMPILER) -c $(dir)mod_interpolate_ctdata.f90

mod_vandermonde.$(obj) : $(dir)mod_vandermonde.f90
	$(COMPILER) -c $(dir)mod_vandermonde.f90

mod_qr.$(obj) : $(dir)mod_qr.f90
	$(COMPILER) -c $(dir)mod_qr.f90

mod_gauss_filter.$(obj) : $(dir)mod_gauss_filter.f90
	$(COMPILER) -c $(dir)mod_gauss_filter.f90

filter.$(obj) : $(dir)filter.f90
	$(COMPILER) -c $(dir)filter.f90

mod_median.$(obj) : $(dir)mod_median.f90
	$(COMPILER) -c $(dir)mod_median.f90

image_process.$(obj) : $(dir)image_process.f90
	$(COMPILER) -c $(dir)image_process.f90

fftpack_module.$(obj) : $(dir)/FFT/source_files/fftpack_module.f90
	$(COMPILER) -c $(dir)/FFT/source_files/fftpack_module.f90

fast_fourier_transform_module.$(obj) : $(dir)/FFT/source_files/fast_fourier_transform_module.f90
	$(COMPILER) -c $(dir)/FFT/source_files/fast_fourier_transform_module.f90

clean:
	rm *.o  *.mod
