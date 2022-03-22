@ Author : Sukrit Jaiswal

import numpy as np
from sympl import Stepper, get_constant
from sympl import initialize_numpy_arrays_with_properties
from datetime import timedelta
from array import array
timestep = timedelta(hours=1)

class surfaceflux(Stepper):
 
    input_properties = {
                
                'air_temperature': {
                'dims': ['mid_levels', '*'],
                'units': 'degK ',
                },

                'specific_humidity': {
                'dims': ['mid_levels', '*'],
                'units': 'kg/kg',
                },

                'air_pressure': {
                'dims': ['mid_levels', '*'],
                'units': 'Pa',
                },

                'northward_wind': {
                'dims': ['mid_levels', '*'],
                'units': 'm s^-1',
                },

                'eastward_wind': {
                'dims': ['mid_levels', '*'],
                'units': 'm s^-1',
                },

                'surface_air_pressure': {
                'dims': ['*'],
                'units': 'Pa',
                },

                'surface_temperature': {
                'dims': ['*'],
                'units': 'degK',
                },

                'surface_specific_humidity': {
                'dims': ['*'],
                'units': 'kg/kg',
                },
                
                'air_pressure_on_interface_levels': {
                'dims': ['interface_levels', '*'],
                'units': 'Pa',
                },
    }
  
 
    diagnostic_properties = {
          
                'richardson_number': {
                'dims': ['*'],'units': 'dimensionless',
                },
 
                'virtual_potential_temperature_at_lowest_level': {
                'dims': ['*'],'units': 'degK'  
                },

                'virtual_potential_temperature_at_roughness_length': {
                'dims': ['*'],'units': 'degK'  
                },
 
                'potential_temperature_at_lowest_level': {
                'dims': ['*'],'units': 'degK'  
                },
 
                'potential_temperature_at_roughness_length': {
                'dims': ['*'],'units': 'degK'  
                },
 
                'surface_upward_sensible_heat_flux': {
                'dims': ['*'],'units': 'W m^-2'
                },

                'surface_upward_latent_heat_flux': {
                'dims': ['*'],'units': 'm^3'
                },

                'height_of_lowest_level':{
                'dims': ['*'],'units': 'm'    
                },

                'drag_coefficient': {
                'dims': ['*'],'units': 'NaN'    
                }, 
    }
 
    output_properties = {
          
                'surface_stress': {
                'dims': ['*'],'units': 'N/m^2'
                },
        
                'surface_north_wind_stress': {
                'dims': ['*'],'units': 'N/m^2'
                },
          
                'surface_east_wind_stress': {
                'dims': ['*'],'units': 'N/m^2'
                },

                'air_temperature': {
                'dims': ['mid_levels', '*'], 'units': 'degK ',
                },

                'specific_humidity': {
                'dims': ['mid_levels', '*'], 'units': 'kg/kg',
                },

                'north_wind_speed': {
                'dims': ['*'], 'units': 'm s^-1',
                },

                'east_wind_speed': {
                'dims': ['*'], 'units': 'm s^-1',
                },

                'new_wind_speed': {
                'dims': ['*'], 'units': 'm s^-1',
                },

                'northward_wind': {
                'dims': ['mid_levels', '*'],
                'units': 'm s^-1',
                },

                'eastward_wind': {
                'dims': ['mid_levels', '*'],
                'units': 'm s^-1',
                },

    }
 
    def _update_constants(self):
        
          """
          Args:
          roughness_length:
            The height at which the wind speed theoretically becomes zero 
            in the absence of wind-slowing obstacles under neutral conditions
          Karman_constant:
            A constant value that is used in the drag coefficient calculation.
          viscosity_coefficient:
            The drag coefficient that is used to calculate
            richardson number
          """
        
          self._g = get_constant('gravitational_acceleration', 'm s^-2')
          self._l = get_constant('latent_heat_of_vaporization_of_water', 'J kg^-1')
          self._cp = get_constant('heat_capacity_of_dry_air_at_constant_pressure', 'J kg^-1 K^-1')
          self._c = 3.21*10**(-5)   # roughness length
          self._Ric = 1
          self._kappa = 0.4   # von Karman constant
          self._v = 10**(-2)    # viscosity_coefficient
          self._ps = 1000   # surface_density
          self._pa = 1001   # air_density
          self._Rd = 287.04
          self._P0 = 100000          
          self._L = 2260 #kJ/kg

 

    def __init__(self,  **kwargs):
    
        self._update_constants()
        super(surfaceflux, self).__init__(**kwargs)
 
 
    def array_call(self, raw_state, timestep):
        
        '''
        Calculates sensible flux and latent heat flux and evaporation returns
        them after timestep.
        '''

        num_cols = raw_state['air_temperature'].shape[1]

        new_state = initialize_numpy_arrays_with_properties(
            self.output_properties, raw_state, self.input_properties
        )

        diagnostics = initialize_numpy_arrays_with_properties(
            self.diagnostic_properties, raw_state, self.input_properties
        )

        new_state['air_temperature'][:] = raw_state['air_temperature']
        new_state['specific_humidity'][:] = raw_state['specific_humidity']
        new_state['northward_wind'][:] = raw_state['northward_wind']
        new_state['eastward_wind'][:] = raw_state['eastward_wind']

        for col in range(num_cols):
            
            air_temperature = raw_state['air_temperature'][:,col]
            surface_temperature = raw_state['surface_temperature'][col]
            air_pressure = raw_state['air_pressure'][:,col]
            surface_pressure = raw_state['surface_air_pressure'][col]
            air_pressure_int = raw_state['air_pressure_on_interface_levels'][:,col]
            specific_humidity = raw_state['specific_humidity'][:,col]
            surface_specific_humidity = raw_state['surface_specific_humidity'][col]
            north_wind = raw_state['northward_wind'][:,col]
            east_wind = raw_state['eastward_wind'][:,col]

            north_wind_int = 0.5*(north_wind[1:] + north_wind[:-1])
            east_wind_int = 0.5*(east_wind[1:] + east_wind[:-1])

            wind_int = np.sqrt(np.power(north_wind_int,2) + np.power(east_wind_int,2))

            air_temperature_int = 0.5*(air_temperature[1:] + air_temperature[:-1])

            specific_humidity_int = 0.5*(specific_humidity[1:] + specific_humidity[:-1])
            
            rho = air_pressure_int[1:-1]/(self._Rd*(1 + 0.608*specific_humidity_int)*air_temperature_int)

            print(rho[0])

            north_wind_speed = north_wind[0]
            east_wind_speed = east_wind[0]

            n = len(air_temperature_int)
        
            pot_virt_temp = air_temperature[0]*np.power((1 - 3.21*10**(-5)),self._Rd/self._cp)*(1 + 0.61*specific_humidity[0])

            pot_virt_temp_surf = surface_temperature*(1 + 0.61*surface_specific_humidity)

#            pot_virt_temp_surf = surface_temperature*np.power((100000/surface_pressure),self._Rd/self._cp)*(1 + 0.61*surface_specific_humidity)

            z = np.zeros(n)
#            z[0] = (self._Rd*(1 + 0.61*specific_humidity[0])*air_temperature[0]/self._g)*np.log(air_pressure_int[0]/air_pressure[0])   #Hydrostatic equation

            z[0] = (surface_pressure - air_pressure[0])/(rho[0]*self._g)

            #x = np.zeros(n)
            #x[0] = (surface_pressure - air_pressure[0])/(rho[0]*self._Rd)
#           height from gas equation
#           print(x[0])

            height_of_lowest_level = z[0]    

            print(z[0])

            potential_temperature_at_lowest_level = surface_temperature*(100000/air_pressure[0])     #Pressure ratio
            potential_temperature_at_roughness_length = air_temperature[0]*(100000/surface_pressure)**(0.286)    #(R/cp = 0.286 for air)        

            #virtual_potential_temperature_at_lowest_level = potential_temperature_at_lowest_level*(1 + 0.61*specific_humidity[0])
            #virtual_potential_temperature_at_roughness_length = potential_temperature_at_roughness_length*(1 + 0.61*specific_humidity[0])
        
 
            virtual_potential_temperature_at_lowest_level = pot_virt_temp
            virtual_potential_temperature_at_roughness_length = pot_virt_temp_surf

            richardson_number = self._g*z[0]*(virtual_potential_temperature_at_lowest_level - virtual_potential_temperature_at_roughness_length)/(virtual_potential_temperature_at_roughness_length*wind_int[0]*wind_int[0]) #richardson number
            # height_of_lowest_level
 
            if richardson_number < 0 :
                drag_coefficient = (self._kappa)**2*(np.log(height_of_lowest_level/self._c))**(-2)          
            elif 0 < richardson_number < 1 :
                drag_coefficient = (self._kappa)**2*(np.log(height_of_lowest_level/self._c))**(-2)*(1-richardson_number)**2          
            else :
                drag_coefficient = 0
 
  
            wind_speed = np.sqrt(pow(north_wind[0], 2) + pow(east_wind[0], 2))
        
            m = (surface_pressure - air_pressure[0])/self._g
            SHF0 = m*self._cp*air_temperature[0]
            LHF0 = m*self._L*specific_humidity[0]
            Q0 = SHF0 + LHF0

            surface_stress = rho[0]*wind_speed*wind_speed*drag_coefficient

            surface_north_wind_stress = rho[0]*north_wind_speed*wind_speed*drag_coefficient
            surface_east_wind_stress = rho[0]*east_wind_speed*wind_speed*drag_coefficient

            surface_upward_sensible_heat_flux = rho[0]*self._cp*wind_speed*drag_coefficient*(surface_temperature - air_temperature[0])*height_of_lowest_level
 
            surface_upward_latent_heat_flux = self._L*rho[0]*wind_speed*drag_coefficient*(surface_specific_humidity - specific_humidity[0])*height_of_lowest_level #evaporation
 
            new_state['air_temperature'][0,col] = raw_state['air_temperature'][0,col] + timestep.total_seconds()*surface_upward_sensible_heat_flux/(height_of_lowest_level*rho[0]*self._cp)  #temperature change

#            new_state['air_temperature'][0,col] = raw_state['air_temperature'][0,col]/(1 + timestep.total_seconds()*surface_upward_sensible_heat_flux/(height_of_lowest_level*rho[0]*self._Rd))  #temperature change


#            print(surface_upward_sensible_heat_flux)
#            print(surface_upward_latent_heat_flux)            
#            print(timestep.total_seconds())

            new_state['specific_humidity'][0,col] = raw_state['specific_humidity'][0,col] + timestep.total_seconds()*surface_upward_latent_heat_flux/(height_of_lowest_level*self._cp*rho[0])
        
#            new_state['specific_humidity'][0,col] = raw_state['specific_humidity'][0,col]/(1 + timestep.total_seconds()*surface_upward_latent_heat_flux/(height_of_lowest_level*self._Rd*rho[0]))


            print(east_wind_speed)
            east_wind_speed = east_wind_speed/(1 + timestep.total_seconds()*surface_east_wind_stress/(height_of_lowest_level*rho[0]))

            north_wind_speed = north_wind_speed/(1 + timestep.total_seconds()*surface_north_wind_stress/(height_of_lowest_level*rho[0]))

#            print(surface_east_wind_stress)
#            print(rho[0])
#            print(height_of_lowest_level)
            
            
# Validation of code

            m = (surface_pressure - air_pressure[0])/self._g
            SHF = m*self._cp*new_state['air_temperature'][0,col]
            LHF = m*self._L*new_state['specific_humidity'][0,col]
            Q = SHF + LHF

            print(m)
            print(Q - Q0)
            print((Q - Q0)/timestep.total_seconds())

            print(surface_upward_sensible_heat_flux)
            print(surface_upward_latent_heat_flux)
            print(surface_upward_sensible_heat_flux + surface_upward_latent_heat_flux)

            new_wind_speed = np.sqrt(pow(north_wind_speed, 2) +
                              pow(east_wind_speed, 2))

            north_wind[0] = north_wind_speed  # new speed
            east_wind[0] = east_wind_speed

            #new_state['air_temperature'][:,col] = air_temperature
            #new_state['specific_humidity'][:,col] = specific_humidity
            new_state['northward_wind'][:,col] = north_wind
            new_state['eastward_wind'][:,col] = east_wind
            

 
            diagnostics = {

                'richardson_number': richardson_number,
  
                'virtual_potential_temperature_at_lowest_level': virtual_potential_temperature_at_lowest_level,

                'virtual_potential_temperature_at_roughness_length': virtual_potential_temperature_at_roughness_length, 
  
                'potential_temperature_at_lowest_level': potential_temperature_at_lowest_level,
 
                'potential_temperature_at_roughness_length': potential_temperature_at_roughness_length,

                'surface_upward_sensible_heat_flux': surface_upward_sensible_heat_flux,

                'surface_upward_latent_heat_flux': surface_upward_latent_heat_flux,

                'height_of_lowest_level': height_of_lowest_level,

                'drag_coefficient': drag_coefficient,
            }
 
        return diagnostics, new_state
 
 
from sympl import PlotFunctionMonitor
from climt import SimplePhysics, get_default_state
#from datetime import timedelta
#timestep = timedelta(hours=1)


def plot_function(fig, state):
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(
        state['specific_humidity'].values.flatten(),
        state['air_pressure'].values.flatten(), '-o')
    ax.axes.invert_yaxis()
    ax.set_ylim(1e5, 100.)
    ax.set_xlabel('kg/kg')
    ax.set_ylabel('Pa')
    ax.set_title('Specific Humidity')


monitor = PlotFunctionMonitor(plot_function)
surface_flux = surfaceflux()
state = get_default_state([surface_flux])
state['eastward_wind'].values[:] = 10.
state['northward_wind'].values[:] = 10.
state['specific_humidity'].values[:] = 0.012
state['surface_specific_humidity'].values[:] = 0.0125
state['air_pressure'].values[:] = 101204.93
state['air_pressure_on_interface_levels'].values[:] = 101325.93
state['surface_air_pressure'].values[:] = 101325.0
state['air_temperature'].values[:] = 288.0    #air temperature decreases by 6.5 degrees C for every 1000 meters you gain.
state['surface_temperature'].values[:] = 288.065
timestep = timedelta(hours=1)
 
for i in range(1):
 
    diagnostics, new_state = surface_flux(state, timestep)
    print('richardson Number:', diagnostics['richardson_number'].values.item())
    print('surface_north_wind_stress', new_state['surface_north_wind_stress'])
    print('surface_east_wind_stress', new_state['surface_east_wind_stress'])
    print('Sensible Heat Flux:', diagnostics['surface_upward_sensible_heat_flux'])
    print('Latent Heat Flux:', diagnostics['surface_upward_latent_heat_flux'])
    print('drag_coefficient:', diagnostics['drag_coefficient'])
    print('New air temperature:', new_state['air_temperature'])
    print('New specific humidity:', new_state['specific_humidity'])
    print('New wind speed:', new_state['new_wind_speed'])
    print('Northward_wind',  new_state['northward_wind'])
    print('Eastward_wind',  new_state['eastward_wind'])

    state.update(diagnostics)
    monitor.store(state)
    state.update(new_state)