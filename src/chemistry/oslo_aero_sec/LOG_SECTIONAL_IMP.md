## Log


### Todo
- FInd out how to edit namelists? 

- Check if I can replace  NPFcoag with NorESM2 variant.
- Can I move aerosect_write2file to aero_model? 
- 


### documentation:
Go through files:

- ~~aeronucl.F90~~
- ~~condtend.F90~~
- ~~aero_sectional.F90~~
- ~~koagsub.F90~~
- ~~aero_model.F90~~ (OBS! ADDED ONE LINE FROM updated noresm (oslo_wetdens(:,:,:) = 0._r8))
- ~~mo_gas_phase_chemdr.F90~~ now identical!! (changes in aero_model.F90 instead in NorESM2) 
- 

### Add option, i.e. make possible to choose between schemes: 
- aeronucl.F90
- condtend.F90
- aero_sectional.F90
- koagsub.F90
- aero_model.F90
- mo_gas_phase_chemdr.F90

### Other challenges:
- Coagulation for nucleating particles --> use alf's code? 
- How to add options? How to reuse as musch code as possible? 
