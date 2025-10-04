import pandas as pd
import requests
import os

def download_sdf(compound_name, output_folder):
    # Construct the URL for the PubChem SDF download
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/ {compound_name}/record/SDF?record_type=3d"
    
    # Make the HTTP request
    response = requests.get(url)
    
    # Check if the request was successful
    if response.status_code == 200:
        # Define the output file path
        output_file = os.path.join(output_folder, f"{compound_name}.sdf")
        
        # Write the content to a file
        with open(output_file, 'wb') as file:
            file.write(response.content)
        print(f"Downloaded {compound_name}.sdf")
    else:
        print(f"Failed to download {compound_name}.sdf")

def main():
    # Define the path to the Excel file, the output folder, sheet name, and column name
    excel_file = 'extra clean.xlsx'     # Replace with your Excel file path
    output_folder = 'sdf_files'         # Replace with your desired output folder
    sheet_name = 'UP makalah'           # Replace with your sheet name
    column_name = 'Metabolite name'     # Replace with your column name
    
    # Create the output folder if it doesn't exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    # Read the specified sheet from the Excel file
    try:
        df = pd.read_excel(excel_file, sheet_name=sheet_name)
    except ValueError as e:
        print(f"Error: {e}")
        print(f"Please check if the sheet '{sheet_name}' exists in the Excel file.")
        return
    
    # Print the column names to check
    print("Column names in the Excel file:", df.columns)
    
    # Check if the required column exists
    if column_name not in df.columns:
        print(f"Error: The column '{column_name}' does not exist in the specified sheet.")
        return
    
    # Iterate over the specified column
    for compound_name in df[column_name]:
        download_sdf(compound_name, output_folder)

if __name__ == "__main__":
    main()
