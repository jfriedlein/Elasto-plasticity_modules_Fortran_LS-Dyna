! -----------MODULE ----------------------------

      module module_dsyevj3

c Usage
       interface spectral_decomp
           module procedure dsyevj3
       end interface
      
      contains

!      ------BEGIN FUNCTIONS-------------------------------------
       include './dsyevj3.f'
c
      end module module_dsyevj3
