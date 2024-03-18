#%%
from aimsgb import GrainBoundary, Grain

#%%
s_input = Grain.from_mp_id("mp-5229")
#%%
gb = GrainBoundary([1, 1, 0], 3, [1, -1, 2], s_input)
structure = Grain.stack_grains(gb.grain_a, gb.grain_b, direction=gb.direction)
# %%
