import pandas as pd
from CoreComputation.EngineClass import engine
from CoreComputation.PropfeedClass import propfeed
from CoreComputation.VehicleClass import vehicle
from rocketcea.cea_obj import CEA_Obj
import numpy as np
import time

start = time.time()

def evaluate(config: dict,outputs: list) -> dict:
    cea_results = CEA_Obj(oxName=config["oxName"], fuelName=config["fuelName"])
    Engine = engine(ceaObj=cea_results, OF=config["OF"], Pc_atm=config["Pc"], M_dot=config["mdot"], Thrust=config["Thrust"])
    PropFeed = propfeed(engine=Engine, burnTime=config["burnTime"], vehicleRadius=config["vehicleRadius"])
    Vehicle = vehicle(engine=Engine, propfeed=PropFeed, vehicleDryMass=config["vehicleDryMass"],
                      vehicleSectionLength=config["vehicleSectionLength"], NoseconeHeight=config["NoseconeHeight"],
                      finLength=config["finLength"], totalFinArea=config["totalFinArea"])

    output_map = {
        "Isp": lambda: Engine.Isp,
        "Thrust": lambda: Engine.Thrust,
        "Mdot": lambda: Engine.M_dot,
        "Fuel Volume": lambda: PropFeed.keroVolume,
        "Oxidizer Volume": lambda: PropFeed.LOXVolume,
        "Launch TWR": lambda: Vehicle.FlightSimulation['LaunchTWR'],
        "Apogee": lambda: Vehicle.FlightSimulation['Apogee']
        # fill out all the properties
    }

    return {key: output_map[key]() for key in outputs}


def twoDimensionalSweep(initialConfig, varOne, rangeOne, varTwo, rangeTwo, outputList):
    rows = []
    for v1 in rangeOne:
        for v2 in rangeTwo:
            config = {**initialConfig, varOne: v1, varTwo: v2}
            result = evaluate(config, outputList)
            row = {varOne: v1, varTwo: v2}
            row.update(result)
            rows.append(row)
    return pd.DataFrame(rows)


config = {
    "OF": 5,
    "Pc": 20,
    "mdot": 2,
    "Thrust": None,
    "burnTime": 20,
    "vehicleRadius": 2,
    "oxName": 'LOX',
    "fuelName":'RP_1',
    "vehicleDryMass": 55,
    "vehicleSectionLength": 10,
    "NoseconeHeight": 1,
    "finLength": 1,
    "totalFinArea": 1
}

outputs = ["Isp", "Thrust", "Oxidizer Volume", "Launch TWR", "Apogee"]
of_values = np.linspace(2.0, 3.0, 20)
pc_values = np.linspace(200, 400, 20)

print(evaluate(config,outputs))
# results = twoDimensionalSweep(config, "OF", of_values, "Pc", pc_values, outputs)
# print(results)

end = time.time()
print(f"Evaluation runtime {end - start} seconds")
