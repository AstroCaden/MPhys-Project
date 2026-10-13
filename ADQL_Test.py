import pyvo as vo
import csv
import pathlib
import tarfile

from astropy.utils.data import download_file
import matplotlib.pyplot as plt
import shutil
from astroquery.eso import Eso
import pandas as pd
import requests
from astroquery.gaia import Gaia
from astropy.io import fits
import os
import numpy as np
from scipy.signal import find_peaks


https://ftp.eso.org/pub/dfs/pipelines/instruments/uves/uves-kit-6.5.3-6.tar.gz
tar xzf uves-kit-6.5.3-6.tar.gz
cd uves-kit-6.5.3-6
./install_pipeline


def downloading_data():
    
    cone_radius = 0.1

    gaia_service = vo.dal.TAPService(
        "https://gea.esac.esa.int/tap-server/tap"
    )

    adql_query = """
    select
        bm.source_id,
        gs.ra,
        gs.dec,
        bm.m1,
        bm.m1_lower,
        bm.m1_upper,
        bm.m2,
        bm.m2_lower,
        bm.m2_upper,
        nss.nss_solution_type,
        nss.significance
    from gaiadr3.binary_masses as bm
    join gaiadr3.nss_two_body_orbit as nss
        on bm.source_id = nss.source_id
    join gaiadr3.gaia_source as gs
        on bm.source_id = gs.source_id
    where bm.m2_lower > bm.m1_upper
      and bm.m2_lower > 1.5
      and bm.m2_lower < 3
      and nss.nss_solution_type = 'SB1'
      and nss.significance > 20
    """

    query_results = gaia_service.search(adql_query)
    results_table = query_results.to_table()

    print(results_table)
    print(f"Number of candidates: {len(results_table)}")


    with open('Gaia_Data_DR3.csv', 'w', newline='') as f:
        writer = csv.writer(f)

        writer.writerow(results_table.colnames)
        writer.writerows(results_table.iterrows())

    eso = Eso()

    eso_service = vo.dal.TAPService("https://archive.eso.org/tap_obs")

    uves_table = []

    for target in results_table:
        ra = target["ra"]
        dec = target["dec"]
        source_id = target["source_id"]

        adql_query = f"""
        select *
        from ivoa.ObsCore
        where instrument_name = 'UVES'
          and intersects(
              circle('ICRS', {ra}, {dec}, {cone_radius}),
              s_region
          ) = 1
        """

        eso_query = eso_service.search(adql_query).to_table()

        print("UVES", source_id, len(eso_query))

        if len(eso_query) > 0:
            eso_query["GAIA_SOURCE_ID"] = source_id
            uves_table.append(eso_query)

    with open('UVES_data.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(uves_table[0].colnames)

        for table in uves_table:
            writer.writerows(table)


    feros_table = []

    for target in results_table:
        ra = target["ra"]
        dec = target["dec"]
        source_id = target["source_id"]

        adql_query = f"""
        select *
        from ivoa.ObsCore
        where instrument_name = 'FEROS'
          and intersects(
              circle('ICRS', {ra}, {dec}, {cone_radius}),
              s_region
          ) = 1
        """

        eso_query = eso_service.search(adql_query).to_table()

        print("FEROS", source_id, len(eso_query))

        if len(eso_query) > 0:
            eso_query["GAIA_SOURCE_ID"] = source_id
            feros_table.append(eso_query)

    with open('FEROS_data.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(feros_table[0].colnames)

        for table in feros_table:
            writer.writerows(table)


    datatype = "both"

    if datatype == "feros":


        os.makedirs("feros_data", exist_ok=True)

        for target in feros_table:
            target_id = str(target["GAIA_SOURCE_ID"][0])
            os.makedirs(f"feros_data/{target_id}", exist_ok=True)
            for dp_id in target["dp_id"]:
                filename = str(dp_id).replace(":", "-") + ".fits"
                url = f"https://dataportal.eso.org/dataPortal/file/{dp_id}"

                response = requests.get(url)
                response.raise_for_status()

                with open(pathlib.Path(f"feros_data/{target_id}") / filename, "wb") as f:
                    f.write(response.content)

    elif datatype == "uves":

        os.makedirs("uves_data", exist_ok=True)

        for target in uves_table:
            target_id = str(target["GAIA_SOURCE_ID"][0])
            os.makedirs(f"uves_data/{target_id}", exist_ok=True)
            for dp_id in target["dp_id"]:
                filename = str(dp_id).replace(":", "-") + ".fits"
                url = f"https://dataportal.eso.org/dataPortal/file/{dp_id}"

                response = requests.get(url)
                response.raise_for_status()

                with open(pathlib.Path(f"uves_data/{target_id}") / filename, "wb") as f:
                    f.write(response.content)



    elif datatype == "both":
        os.makedirs("uves_data", exist_ok=True)

        for target in uves_table:
            target_id = str(target["GAIA_SOURCE_ID"][0])
            os.makedirs(f"uves_data/{target_id}", exist_ok=True)
            for dp_id in target["dp_id"]:
                filename = str(dp_id).replace(":", "-") + ".fits"
                url = f"https://dataportal.eso.org/dataPortal/file/{dp_id}"

                response = requests.get(url)
                response.raise_for_status()

                with open(pathlib.Path(f"uves_data/{target_id}") / filename, "wb") as f:
                    f.write(response.content)



        os.makedirs("feros_data", exist_ok=True)

        for target in feros_table:
            target_id = str(target["GAIA_SOURCE_ID"][0])
            os.makedirs(f"feros_data/{target_id}", exist_ok=True)
            for dp_id in target["dp_id"]:
                filename = str(dp_id).replace(":", "-") + ".fits"
                url = f"https://dataportal.eso.org/dataPortal/file/{dp_id}"

                response = requests.get(url)
                response.raise_for_status()

                with open(pathlib.Path(f"feros_data/{target_id}") / filename, "wb") as f:
                    f.write(response.content)
        



def barycentric_corrections():

    uves_array = os.listdir("uves_data")
    for folder in uves_array:
        folder_data = os.listdir(f"uves_data/{folder}")
        for file in folder_data:
            hdul = fits.open(f"uves_data/{folder}/{file}")
            header = hdul[0].header
            
            if header["SPECSYS"] == "TOPOCENT":
                print(f"{file} is topocentric, performing barycentric corrections")
                






                
            elif header["SPECSYS"] == "BARYCENT":
                print(f"{file} is barycentric, skipping")
                shutil.copy(f"{file}.fits", f"{file}_bary.fits")


                
            else:
                print("SYSSPEC not recognised, skipping")
                continue
        


    feros_array = os.listdir("feros_data")

    for folder in feros_array:
        folder_data = os.listdir(f"feros_data/{folder}")

        for file in folder_data:
            file_path = f"feros_data/{folder}/{file}"

            with fits.open(file_path) as hdul:
                header = hdul[0].header
                specsyst = header.get("SPECSYS")

            if specsyst == "TOPOCENT":
                print(f"{file} is topocentric, performing barycentric corrections")

            elif specsyst == "BARYCENT":
                print(f"{file} is barycentric, copying file")

                name, ext = os.path.splitext(file)
                new_path = f"feros_data/{folder}/{name}_bary{ext}"

                shutil.copy(file_path, new_path)

            else:
                print(f"{file}: SPECSYS not recognised, skipping")
                continue
            




def plot_spectra():
    
    uves_array = os.listdir("uves_data")
    os.makedirs("spectra_UVES", exist_ok=True)

    os.makedirs("combined_spectra", exist_ok=True)
    os.makedirs("combined_spectra/FEROS", exist_ok=True)
    os.makedirs("combined_spectra/UVES", exist_ok=True)

    feros_array = os.listdir("feros_data")
    os.makedirs("spectra_FEROS", exist_ok=True)

    grid_spectra_uves = []
    grid_spectra_feros = []

    if len(uves_array) == 0:
        print("No UVES data found. Skipping")

    else:

    
        print(uves_array)
        for folder in uves_array:
            
            folder_data = os.listdir(f"uves_data/{folder}")
            os.makedirs(f"spectra_UVES/{folder}", exist_ok=True)
            spectra = []
            wavelengths = []
            for file in folder_data:
                
                hdul = fits.open(f"uves_data/{folder}/{file}")
                data = hdul[1].data

                wavelength = data["WAVE"][0]
                flux = data["FLUX"][0]

                
                plt.plot(data["WAVE"][0], data["FLUX"][0], linewidth=1)
                plt.xlabel("Wavelength (Å)")
                plt.ylabel("Flux (10$^{-16}$ erg cm$^{-2}$ s$^{-1}$ Å$^{-1}$)")
                #plt.xlim(6000, 8000)
                plt.savefig(f"spectra_UVES/{folder}/{file}_spectra_UVES.png")
                plt.close()
                print(f"UVES: {file} spectra saved.")
                wavelengths.append(data["WAVE"][0])
                spectra.append(data["FLUX"][0])

            
            common_wave = wavelengths[0]
            spectra = [np.interp(common_wave, wave, flux) for wave, flux in zip(wavelengths, spectra)]
            average_flux = np.mean(spectra, axis=0)

            grid_spectra_uves.append((folder, common_wave, average_flux))

            plt.plot(common_wave, average_flux, linewidth=1)
            plt.xlabel("Wavelength (Å)")
            plt.ylabel("Flux (10$^{-16}$ erg cm$^{-2}$ s$^{-1}$ Å$^{-1}$)")
            plt.savefig(f"combined_spectra/UVES/average_spectrum_{folder}.png")
            plt.close()
            print(f"UVES: Combined spectra {folder} saved.")



    if len(feros_array) == 0:
        print("No feros data found. Skipping")

    else:

    
        print(feros_array)
        for folder in feros_array:
            
            folder_data = os.listdir(f"feros_data/{folder}")
            os.makedirs(f"spectra_FEROS/{folder}", exist_ok=True)
            spectra = []
            wavelengths = []
            for file in folder_data:
                
                hdul = fits.open(f"feros_data/{folder}/{file}")
                data = hdul[1].data


                
                plt.plot(data["WAVE"][0], data["FLUX"][0], linewidth=1)
                plt.xlabel("Wavelength (Å)")
                plt.ylabel("Flux (ADU)")
                plt.ylim(np.nanpercentile(data["FLUX"][0], [1, 99]))
                #plt.xlim(6000, 8000)
                plt.savefig(f"spectra_FEROS/{folder}/{file}_spectra_FEROS.png")
                plt.close()
                print(f"FEROS: {file} spectra saved.")
                wavelengths.append(data["WAVE"][0])
                spectra.append(data["FLUX"][0])

            
            common_wave = wavelengths[0]
            spectra = [np.interp(common_wave, wave, flux) for wave, flux in zip(wavelengths, spectra)]
            average_flux = np.mean(spectra, axis=0)

            grid_spectra_feros.append((folder, common_wave, average_flux))

            plt.plot(common_wave, average_flux, linewidth=1)
            plt.xlabel("Wavelength (Å)")
            plt.ylabel("Flux (ADU)")
            plt.ylim(np.nanpercentile(average_flux, [1, 99]))
            plt.savefig(f"combined_spectra/FEROS/average_spectrum_{folder}.png")
            plt.close()
            print(f"FEROS: Combined spectra {folder} saved.")





    flux_labels = {
        "UVES": "Flux (10$^{-16}$ erg cm$^{-2}$ s$^{-1}$ Å$^{-1}$)",
        "FEROS": "Flux (ADU)"
    }

    wave_ranges = {
        "UVES": [("blue", 0, 5760), ("red", 5760, np.inf)],
        "FEROS": [("full", 0, np.inf)]
    }

    for instrument, grid_spectra in [
        ("UVES", grid_spectra_uves),
        ("FEROS", grid_spectra_feros)
    ]:

        for range_name, wave_lo, wave_hi in wave_ranges[instrument]:

            masks = [(wave >= wave_lo) & (wave < wave_hi) & (flux != 0) for _, wave, flux in grid_spectra]

            x_min = min(np.min(wave[mask]) for (_, wave, flux), mask in zip(grid_spectra, masks))
            x_max = max(np.max(wave[mask]) for (_, wave, flux), mask in zip(grid_spectra, masks))
            y_min = min(np.nanpercentile(flux[mask], 1) for (_, wave, flux), mask in zip(grid_spectra, masks))
            y_max = max(np.nanpercentile(flux[mask], 99) for (_, wave, flux), mask in zip(grid_spectra, masks))

            for start in range(0, len(grid_spectra), 8):

                group = grid_spectra[start:start + 8]
                nrows = (len(group) + 1) // 2

                fig, axes = plt.subplots(nrows, 2, figsize=(14, 3 * nrows), sharex=True, sharey=True, squeeze=False)
                axes = axes.flatten()

                for ax, (folder, wavelength, flux) in zip(axes, group):
                    ax.plot(wavelength, flux, linewidth=0.7)
                    ax.set_title(folder, fontsize=9)
                    ax.set_xlim(x_min, x_max)
                    ax.set_ylim(y_min, y_max)
                    ax.grid(True, alpha=0.2)

                for ax in axes[len(group):]:
                    ax.set_visible(False)

                fig.supxlabel("Wavelength (Å)")
                fig.supylabel(flux_labels[instrument])
                fig.tight_layout()

                fig.savefig(
                    f"combined_spectra/{instrument}_spectral_grid_{range_name}_{start // 8 + 1}.png",
                    dpi=300,
                    bbox_inches="tight"
                )
                plt.close(fig)






def candidates_on_CMD():

    
    bg = Gaia.launch_job_async("""
    select top 100000 parallax, phot_g_mean_mag, bp_rp
    from gaiadr3.gaia_source
    where parallax > 6.6667
      and parallax_over_error > 10
      and bp_rp is not null
      and phot_g_mean_mag is not null
    order by random_index
    """).get_results().to_pandas()

    ids = pd.read_csv('Gaia_Data_DR3.csv', dtype={'source_id': 'int64'})['source_id']
    cand = Gaia.launch_job(f"""
    select source_id, parallax, phot_g_mean_mag, bp_rp
    from gaiadr3.gaia_source
    where source_id in ({','.join(map(str, ids))})
    """).get_results().to_pandas()

    for df in (bg, cand):
        df.columns = df.columns.str.lower()
        df['MG'] = df['phot_g_mean_mag'] + 5 * np.log10(df['parallax']) - 10

    plt.figure(figsize=(6, 7))
    plt.scatter(bg['bp_rp'], bg['MG'], s=0.5, c='0.7', lw=0, rasterized=True,
                label='Gaia DR3 sample')
    uves = cand['source_id'].isin(pd.read_csv('UVES_data.csv')['GAIA_SOURCE_ID'])
    feros = cand['source_id'].isin(pd.read_csv('FEROS_data.csv')['GAIA_SOURCE_ID'])

    plt.scatter(cand['bp_rp'][uves & ~feros], cand['MG'][uves & ~feros], s=40, c='#1baf7a', edgecolor='k', zorder=3, label='UVES')
    plt.scatter(cand['bp_rp'][feros & ~uves], cand['MG'][feros & ~uves], s=40, c='#eb6834', edgecolor='k', zorder=3, label='FEROS')
    plt.scatter(cand['bp_rp'][uves & feros], cand['MG'][uves & feros], s=40, c='#2a78d6', edgecolor='k', zorder=3, label='Both')
    plt.gca().invert_yaxis()
    plt.xlabel('BP − RP')
    plt.ylabel('M_G')
    plt.legend()
    plt.savefig("candidates_on_CMD.png")
    plt.show()












#downloading_data()
barycentric_corrections()
#plot_spectra()
#candidates_on_CMD()
