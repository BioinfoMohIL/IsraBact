import time
import json
import pandas as pd
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

# === Setup Selenium Headless Chrome ===
chrome_options = Options()
chrome_options.add_argument("--headless")
driver = webdriver.Chrome(options=chrome_options)

url = "https://virulence.fda.gov/data/dnasequences/?type=p"
driver.get(url)

# === Attente que le tableau charge ===
WebDriverWait(driver, 20).until(
    EC.presence_of_element_located((By.XPATH, '//*[@id="DataTables_Table_0"]/tbody'))
)
time.sleep(2)

data = []

while True:
    rows = driver.find_elements(By.XPATH, '//*[@id="DataTables_Table_0"]/tbody/tr')
    print(f"🧪 Found {len(rows)} rows on this page")

    for row in rows:
        try:
            genome_id = row.find_element(By.XPATH, "./td[1]").text.strip()
            dna_sequence = row.find_element(By.XPATH, "./td[2]").text.strip()
            data.append({"genome_id": genome_id, "dna_sequence": dna_sequence})
        except Exception as e:
            print(f"⚠️ Skipping row due to error: {e}")

    # === Pagination: vérifie si le bouton "Next" est cliquable ===
    try:
        next_button = driver.find_element(By.ID, "DataTables_Table_0_next")
        next_class = next_button.get_attribute("class")
        if "disabled" in next_class:
            print("✅ Reached the last page.")
            break
        else:
            next_button.click()
            time.sleep(1.5)  # laisser le temps de recharger
    except Exception as e:
        print(f"❌ Pagination error: {e}")
        break

driver.quit()

# === Save JSON ===
with open("plasmid_dna_sequences.json", "w") as json_file:
    json.dump(data, json_file, indent=2)

# === Save CSV ===
df = pd.DataFrame(data)
df.to_csv("plasmid_dna_sequences.csv", index=False)

print("✅ Done. Data saved to JSON and CSV.")
