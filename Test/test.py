def state_preserve_query(produced_value: str):
	"""Generate a dictionary that preserves the state of a generic query between simulation steps.
		Input the name of the preserved value as a string.	
	"""
	consumed_query = '''(
				prev!(''' + produced_value + '''),
			)'''
	return {'consumed': consumed_query, 'produced': produced_value, 'function': identity}

def identity(arg):
    return arg


def mass_burn_preserved_properties():
	preserved_properties = ['massBurned', 'burnRate', 'thrustAccel', 'thrustDirection']

	queries = []
	for preserved_property in preserved_properties:
		queries.append(state_preserve_query(preserved_property))
		
	return queries

list1 = [1,2]

list2 = []
list2.extend([3])

list3 = [4,5]

total = list1 + list2 + list3

print(total)

print(mass_burn_preserved_properties())
