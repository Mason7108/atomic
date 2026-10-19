import json
from rdkit import Chem

ptable=Chem.GetPeriodicTable()
order=[('1s',2),('2s',2),('2p',6),('3s',2),('3p',6),('4s',2),('3d',10),('4p',6),('5s',2),('4d',10),('5p',6),('6s',2),('4f',14),('5d',10),('6p',6),('7s',2),('5f',14),('6d',10),('7p',6)]
# Known non-Aufbau ground-state electron configurations (selected exceptions).
exceptions={
 24:{'4s':1,'3d':5},29:{'4s':1,'3d':10},
 41:{'5s':1,'4d':4},42:{'5s':1,'4d':5},44:{'5s':1,'4d':7},45:{'5s':1,'4d':8},46:{'5s':0,'4d':10},47:{'5s':1,'4d':10},
 57:{'4f':0,'5d':1},58:{'4f':1,'5d':1},64:{'4f':7,'5d':1},71:{'4f':14,'5d':1},
 78:{'6s':1,'5d':9},79:{'6s':1,'5d':10},
 89:{'5f':0,'6d':1},90:{'5f':0,'6d':2},91:{'5f':2,'6d':1},92:{'5f':3,'6d':1},93:{'5f':4,'6d':1},96:{'5f':7,'6d':1},
 103:{'6d':0,'7p':1},
}
# For select exceptions, moving electrons from one subshell to another preserves Z.
# For 90-93 orbital filling starts from 7s2 5f; multiple deviations are manually assigned.
full_orbitals=[('1s',2),('2s',2),('2p',6),('3s',2),('3p',6),('4s',2),('3d',10),('4p',6),('5s',2),('4d',10),('5p',6),('6s',2),('4f',14),('5d',10),('6p',6),('7s',2),('5f',14),('6d',10),('7p',6)]

def electronic(z):
 left=z
 vals={}
 for orb,cap in full_orbitals:
  used=min(left,cap); vals[orb]=used; left-=used
  if not left:break
 # Exception definitions assign values for given orbitals, correcting in order via deltas.
 if z in exceptions:
  for orb,val in exceptions[z].items():vals[orb]=val
  # Balance by redistributing from the nearest relevant partially-filled subshell.
  diff=z-sum(vals.values())
  if diff:
   dest={'57':'4f'}
   # For lanthanide and actinide exceptions, our base Aufbau 4f/5f fills early; 5d/6d is changed,
   # and balancing should adjust f first, except Lr where 6d becomes 7p.
   pool='4f' if 57<=z<=71 else '5f' if 89<=z<=103 else '6s' if 76<=z<=79 else '5s' if 41<=z<=47 else '4s'
   # 57 config [Xe]5d1 6s2 means initial 4f1 -> 4f0 & 5d1 (no delta).
   # 90 etc need special full configuration overrides.
   vals[pool]=vals.get(pool,0)+diff
 # Correct special f-block configurations using accurate configurations
 special={
  57:{'6s':2,'4f':0,'5d':1},58:{'6s':2,'4f':1,'5d':1},64:{'6s':2,'4f':7,'5d':1},71:{'6s':2,'4f':14,'5d':1},
  89:{'7s':2,'5f':0,'6d':1},90:{'7s':2,'5f':0,'6d':2},91:{'7s':2,'5f':2,'6d':1},92:{'7s':2,'5f':3,'6d':1},93:{'7s':2,'5f':4,'6d':1},96:{'7s':2,'5f':7,'6d':1},
  103:{'7s':2,'5f':14,'6d':0,'7p':1}
 }
 if z in special:
  vals.update(special[z])
 shells=[0]*7
 for orb,v in vals.items():shells[int(orb[0])-1]+=v
 assert sum(shells)==z, (z,vals,shells)
 while shells and shells[-1]==0:shells.pop()
 return shells,' '.join(f'{orb}{vals[orb]}' for orb,_ in full_orbitals if vals.get(orb,0)>0)

alkali={3,11,19,37,55,87}
alkaline={4,12,20,38,56,88}
halogens={9,17,35,53,85,117}
nobles={2,10,18,36,54,86,118}
metalloids={5,14,32,33,51,52}
nonmetals={1,6,7,8,15,16,34}
post={13,31,49,50,81,82,83,84,113,114,115,116}
gases={1,2,7,8,9,10,17,18,36,54,86}
liquids={35,80}
# Isotopes to explore (observed known nuclides).
isotopes={1:[1,2,3],2:[3,4],3:[6,7],4:[9,10],5:[10,11],6:[12,13,14],7:[14,15],8:[16,17,18],9:[19],10:[20,21,22],11:[23],12:[24,25,26],13:[27],14:[28,29,30],15:[31,32,33],16:[32,33,34,36],17:[35,37],18:[36,38,40],19:[39,40,41],20:[40,42,43,44,46,48],26:[54,56,57,58],29:[63,65],47:[107,109],79:[197],92:[234,235,238]}

def category(z):
 if z in alkali:return 'Alkali metal'
 if z in alkaline:return 'Alkaline earth metal'
 if z in halogens:return 'Halogen'
 if z in nobles:return 'Noble gas'
 if 57<=z<=71:return 'Lanthanide'
 if 89<=z<=103:return 'Actinide'
 if z in metalloids:return 'Metalloid'
 if z in nonmetals:return 'Nonmetal'
 if z in post:return 'Post-transition metal'
 if (21<=z<=30) or (39<=z<=48) or (72<=z<=80) or (104<=z<=112):return 'Transition metal'
 return 'Unknown / predicted'

# Approximate group number for common groups; f-block group left undefined.
group_numbers={1:1,2:18}
periods=[2,10,18,36,54,86,118]
def period(z):
 return next(i+1 for i,lim in enumerate(periods) if z<=lim)

def group(z):
 if 57<=z<=71 or 89<=z<=103:return None
 if z==1:return 1
 if z==2:return 18
 # Use conventional long-form main periodic table by row.
 ranges=[(3,4,13,18),(11,12,13,18),(19,30,31,36),(37,48,49,54),(55,56,72,86),(87,88,104,118)]
 for a,b,c,d in ranges:
  if a<=z<=b:return z-a+1
  if c<=z<=d:
   # In periods 2 and 3: p-block 13..18; period 4+ p-block 13..18.
   if a in [3,11]:return z-c+13
   # In periods 6/7 transition metal region is 72..80 (group4..12).
   if a in [55,87] and z<= (80 if a==55 else 112):return 4+(z-c)
   if a in [55,87]:return 13+(z-(81 if a==55 else 113))
   return z-c+13 if z>=31 and a in [19,37] else None
 if 21<=z<=30:return z-18
 if 39<=z<=48:return z-36
 return None

out=[]
for z in range(1,119):
 shell,config=electronic(z)
 try:m=ptable.GetMostCommonIsotope(z)
 except: m=round(ptable.GetAtomicWeight(z))
 if m<z:m=round(ptable.GetAtomicWeight(z))
 out.append(dict(number=z,name=ptable.GetElementName(z),symbol=ptable.GetElementSymbol(z),atomicMass=round(ptable.GetAtomicWeight(z),4),isotope=m,isotopes=isotopes.get(z,[m]),shells=shell,configuration=config,category=category(z),period=period(z),group=group(z),phase=('Gas' if z in gases else 'Liquid' if z in liquids else 'Unknown' if z>109 else 'Solid')))

print('118 elements;',[(x['symbol'],x['shells']) for x in out[:20]])
print('U:',out[91]);print('Og:',out[-1]);
with open('/mnt/data/atom-explorer/elements.json','w') as f:json.dump(out,f,ensure_ascii=False,separators=(',',':'))
