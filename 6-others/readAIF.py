from gemmi import cif

data = cif.read_file("../../database/writing/database/CoREMOF2024DB_public/water/ASR/2016[Cu][tfw]2[ASR]1.aif")
block = data.sole_block()
print(block.name)

item = block.find_pair_item('_units_loading')
print(item.pair)

print("uptake:",list(block.find_loop('_adsorp_amount')))
print("pressure:",list(block.find_loop('_adsorp_pressure')))
print("p0:",list(block.find_loop('_adsorp_p0')))