
RATIO_OF_SPECIFIC_HEATS = 7.0/5.0 # for air

def pressure_after_adiabatic_expansion(starting_pressure,starting_gas_volume,current_gas_volume):
    density_ratio = (starting_gas_volume/current_gas_volume)
    new_pressure = starting_pressure * density_ratio**RATIO_OF_SPECIFIC_HEATS
    return new_pressure

def temperature_after_adiabatic_expansion(starting_temperature,starting_gas_volume,current_gas_volume):
    density_ratio = (starting_gas_volume/current_gas_volume)
    new_temperature = starting_temperature * density_ratio**(1.0 - RATIO_OF_SPECIFIC_HEATS)
    return new_temperature

def water_reaction_force(nozzle_area,velocity_relative_to_nozzle,water_density):
    return water_density*nozzle_area*velocity_relative_to_nozzle**2