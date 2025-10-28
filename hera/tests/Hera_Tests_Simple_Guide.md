# Run Hera JSON Tests — Simple Step‑by‑Step (English)

This guide helps you run the tests on **any computer**.  
Just follow the steps in order. No coding knowledge needed.

---

## What you need before you start
- You have the **Hera project** on your machine.
- You have the **test data folder** (usually named `hera_unittest_data`) somewhere on your disk.
- You know where the **tests folder** is (usually `hera/tests` or sometimes `hera/hera/tests`).

---

## 1) Open a terminal
- On **Linux/macOS**: open the Terminal app.
- On **Windows**: open **PowerShell**.

---

## 2) Tell Hera where the data is (set environment variables)

### Linux / macOS (copy–paste and change the path to your own)
```bash
export HERA_DATA_PATH=/path/to/hera_unittest_data
export HERA_UNITTEST_DATA=/path/to/hera_unittest_data
```

### Windows PowerShell (copy–paste and change the path to your own)
```powershell
$env:HERA_DATA_PATH = "C:/path/to/hera_unittest_data"
$env:HERA_UNITTEST_DATA = "C:/path/to/hera_unittest_data"
```

**Explanation:**  
Hera needs to know where your **test data** lives. Replace `/path/to/hera_unittest_data` (or `C:/path/...`) with the **real folder** on **your** machine.

> Tip: If you are not sure, search your disk for a folder that contains files like `.tif`, `.hgt`, `.shp`, `.parquet` used by the tests.

---

## 3) Go to the tests folder (the folder with the run script)
In the same terminal, change directory to the tests folder.

**Option A (most common):**
```bash
cd /path/to/hera/tests
```

**Option B (if your repo has two `hera` folders):**
```bash
cd /path/to/hera/hera/tests
```

**Explanation:**  
You need to be **inside** the folder that contains the script `run_all_json_tests.sh`.

---

## 4) Run in **Save mode** (creates/updates the expected results)
Use this the **first time**, or when you **want to update** the baseline results.

```bash
PREPARE_EXPECTED_OUTPUT=1 bash run_all_json_tests.sh
```

**What you will see:**  
- The tests run and **save** the expected outputs to disk.
- This becomes your “reference” for future comparisons.

---

## 5) Run in **Compare mode** (checks current results vs expected)
Use this for **normal validation**. It will **compare** current results with what you saved in the previous step.

```bash
PREPARE_EXPECTED_OUTPUT=0 bash run_all_json_tests.sh
```

**What you will see:**  
- The tests run and compare outputs.
- If something changed unexpectedly, you will see **differences** reported.
- If everything matches, you will get a **success** message.

---

## 6) (Optional) One command to do both, one after the other
```bash
bash -lc 'PREPARE_EXPECTED_OUTPUT=1 bash run_all_json_tests.sh && PREPARE_EXPECTED_OUTPUT=0 bash run_all_json_tests.sh'
```

**Explanation:**  
- First it saves the expected outputs, then it immediately compares against them.

---

## 7) Common problems and fixes
- **“No such file or directory”** → Check that you are in the **correct tests folder**.
- **“Data not found”** → Check that `HERA_DATA_PATH` and `HERA_UNITTEST_DATA` are set to the **correct data folder**.
- **Windows**: Make sure you are using **PowerShell** and paths look like `C:/...` (forward slashes are OK).

---

## That’s it!
These are the only steps you need: set the data folder, go to the tests folder, and run **Save** or **Compare** mode.
