import time
import json
import pandas as pd
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

# === Headless Chrome setup ===
chrome_options = Options()
chrome_options.add_argument("--headless")
driver = webdriver.Chrome(options=chrome_options)

# === Page URL ===
url = "https://virulence.fda.gov/data/cds/?type=p"
driver.get(url)

# === Wait for second table to load ===
WebDriverWait(driver, 20).until(
    EC.presence_of_element_located((By.ID, "cd-table"))
)
time.sleep(2)

# === Scrape cd-table ===
data = []
while True:
    rows = driver.find_elements(By.XPATH, '//*[@id="cd-table"]/tbody/tr')
    print(f"🔍 Found {len(rows)} rows in cd-table")

    for row in rows:
        cells = row.find_elements(By.TAG_NAME, "td")
        if len(cells) != 7:
            continue
        entry = {
            "ID": cells[0].text.strip(),
            "Virulence_Gene": cells[1].text.strip(),
            "Locus_Tag": cells[2].text.strip(),
            "Note": cells[3].text.strip(),
            "Product": cells[4].text.strip(),
            "DB_Xref": cells[5].text.strip(),
            "Location": cells[6].text.strip()
        }
        data.append(entry)

    # Pagination (if needed) – this example assumes table fits in one page
    try:
        next_btn = driver.find_element(By.ID, "cd-table_next")
        next_class = next_btn.get_attribute("class")
        if "disabled" in next_class:
            break
        else:
            next_btn.click()
            time.sleep(1.5)
    except:
        break

driver.quit()

# === Save JSON ===
with open("virulence_factors.json", "w") as jf:
    json.dump(data, jf, indent=2)

# === Save CSV ===
df = pd.DataFrame(data)
df.to_csv("virulence_factors.csv", index=False)

print("✅ virulence_factors.json and .csv created successfully.")
