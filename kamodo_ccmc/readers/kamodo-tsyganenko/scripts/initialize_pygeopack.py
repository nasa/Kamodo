

import PyGeopack as gp

def update_geopack_params():
    # Update the PyGeopack parameters
    gp.Params.UpdateParameters()
    print("Parameters updated successfully.")

if __name__ == "__main__":
    update_geopack_params()
