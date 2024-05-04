# Database Manager Tool

## Overview
This Database Manager Tool is designed to facilitate database operations for the LA Crime Map project. It enables administrators to perform various CRUD operations through a command-line interface (CLI), connecting to MySQL databases hosted on AWS RDS. The tool supports operations such as data insertion, updating, deletion, and querying.

## Setup and Installation

### Prerequisites
- Python 3.x
- MySQL connector for Python (`mysql-connector-python`)
- Access to AWS RDS with appropriate credentials

### Installation
1. Ensure Python 3.x is installed on your system.
2. Install the MySQL connector package:


## Configuration
Before running the tool, configure the database connection settings in the script. Update the placeholders in the `db_config` dictionary with the actual database host, port, user, and password.

## Usage
To use the tool, run the script from the command line with the required operation and parameters. Here’s how to perform various operations:

### Single Operations

1. **Insert**:
- Description: Provide all necessary fields in the order defined by your database schema, separated by spaces.
- Example:
  ```
  python db_manager.py insert 123456 "2020-01-08 00:00:00" "2020-01-08 12:30:00" 1230 3 "Southwest" "0377" "624" "BATTERY - SIMPLE ASSAULT" 36 "F" "B" 501 "SINGLE FAMILY DWELLING" "400" "STRONG-ARM" "IC" "Invest Cont" "1100 W 39TH PL" "Main St" 34.0141 -118.2978
  ```

2. **Update**:
- Description: Specify the record identifier followed by the field-value pairs you wish to update, prefixed with `--`.
- Example:
  ```
  python db_manager.py update 123456 --Vict_Age=37 --Status="AA" --Status_Desc="Adult Arrest"
  ```

3. **Delete**:
- Description: Provide the record identifier of the entry you wish to delete.
- Example:
  ```
  python db_manager.py delete 123456
  ```

4. **View**:
- Description: Use this function to view all records or specify criteria to filter the results.
- Example:
  ```
  python db_manager.py view --AREA_NAME="Southwest"
  ```

### Bulk Operations

1. **Bulk Insert**:
- Description: Provide the path to a CSV file containing multiple records that match the database schema.
- Example:
  ```
  python db_manager.py bulk_insert "crime_new.csv"
  ```

2. **Bulk Update**:
- Description: Specify a condition and the field-value pairs to update, reflecting the SQL `UPDATE` syntax.
- Example:
  ```
  python db_manager.py bulk_update "Status='AA', Status_Desc='Adult Arrest' WHERE AREA_NAME='Southwest'"
  ```

3. **Bulk Delete**:
- Description: Provide conditions for the records you wish to delete, reflecting the SQL `DELETE` syntax.
- Example:
  ```
  python db_manager.py bulk_delete "AREA_NAME='Southwest'"
  ```
