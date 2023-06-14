# Omnomnomivores

## Community Dynamics Just aren't a Vibe...

When we use the (what is hopefully) a carbon copy of the original Thompson &
Gonzales model the abundances for species are just weird and seems to be a plant
safe space as they are the only species to persist. See below:

!['normal model'](figures/diagnostics.png)

### Test 1: Set all but plant-plant interactions to zero

Here all interactions are zero-ed except for plant-plant interactions. Here we
can see if the root cause of the sus dynamics is because herbivores and
carnivores are going extinct...

!['only plant-plant'](figures/diagnostics_only_plant-plant.png)

This could be the case... dynamics (shape-wise seem pretty unchanged)

### Test 2: Lets ramp up the strength of interactions

Lets start by setting the scaling parameter to zero.

!['no scaling parameter'](figures/diagnostics_no_scaling.png)

I did this by multiplying the interaction values by (0.33*S) since that should
neutralise the original division term, in theory...

### Test 3: Set all interactions to zero

One last hail mary before we try and throw parameters at the problem. Lets set
_ALL_ interactions to zero and see what happens to the plants.

!['no interactions'](figures/diagnostics_no_interactions.png)

### Test 4: The test that should've been the first test

Double, triple, etc. check that the correct values are used for assigning
interaction strengths...

Looking at the source code for T&G some interesting notes.

- The seem to have an 'upside down trophic structure' in the sense that they
  have 5 'plants', 7 'herbivores', and 11 'carnivores'
- All intraspecific interactions are set to 0.2 (all trophic levels)
- Regarding previous it seems they then 'reset' `carnivore` intraspecific
  interactions to 0.15
- weight is defines as '1/80*3' and 80 = S (since there should be 80 species)
  BUT when constructing the food web based on the numbers given on the first
  food web there are only 23 species...  
- Adding to the previous point nSpecies is set to 30...

What is life???

### Test 5: Lets try scaling of 1/80*3

Since this is the scaling factor used lets give it a try...

!['correct scaling param'](figures/diagnostics_correct_scaling.png)

So the abundance numbers 'match' those from T&G script for what thats worth...

### Test 6: Set all initial abundances for same time

Set all initial abundances for all trophic levels for the first timestamp (no
community assembly)

!['no gradual invasion'](figures/diagnostics_same_init.png)

### Test 7: Invert the trophic level

Lets invert the number of species at the different trophic levels.

!['inverted trophic levels'](figures/diagnostics_inverted_trophic.png)

I mean they persist a bit longer. And also plants are still persisting so it is
for sure the internal growth rate of the non plants at this point?

### Test 8: Inverted and gradual invasion

Lets go with the inverted trophic level but re-introduce gradual invasions.

!['inverted trophic levels' but make it gradual invasions](figures/diagnostics_inverted_gradual.png)

Nope!
