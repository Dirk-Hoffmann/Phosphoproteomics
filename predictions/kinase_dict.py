from urllib import response
import requests as r

url = "https://rest.uniprot.org/uniprotkb/search?query=gene:RPS6KA1"

response = r.post(url)

lines= []

for f in response.iter_lines():



    lines.append(f)

print(lines[0])