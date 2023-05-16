#!/bin/bash

##############################################################################################################
# This script produces a split membrane out of a membrane obtained from CHARMM-GUI website with CHARMM36 FF. #       
# Before using the script, pleasde read the README.txt                                                       #
##############################################################################################################

echo
echo "############################"
echo "### Starting the process ###"
echo "############################"
echo

# Read the final output .gro file of Charmm-gui (step7_10.gro)

echo "### Reading the input file..."

# Get lipid types

residues=`tail -n +3 step7_10.gro | head -n-1 | awk -v FIELDWIDTHS="5 4 99" '{print $2}' | uniq` 
residue_actual=`echo $residues | awk -v lipid_list=$list '
BEGIN{
  lipid_array[1]="POPC"
  lipid_array[2]="POPE"
  lipid_array[3]="POPS"
  lipid_array[4]="POPA"
  lipid_array[5]="PSM"
  lipid_array[6]="CHL1"
  lipid_count=6
}
{
  for(i=1; i<=NF; i++) 
    for(j=1; j<=lipid_count; j++ )
      if($i == lipid_array[j])
        print $i
}'`

# Make index file

    {
echo "keep 0"
for i in $residue_actual
do
  echo -n "r $i |"
done
echo
echo "name 1 MEMB"
echo q
    } | gmx make_ndx -f step7_10.gro -o index.ndx

# Remove the title and box size

awk -F '|' -v OFS='|' '{sub(/^.............../, "& "); print}' step7_10.gro > system.gro

head -n-1 system.gro | tail -n+3 > file.tmp && mv file.tmp system.gro

# Remove water and ions

awk '{	if ( $1 !~ /TIP/ )
		if ( $1 !~ /SOD/ )
		if ( $1 !~ /CLA/ )	   
        printf "%9s %5s %5s %7.3f %7.3f %7.3f\n", $1, $2, $3, $4, $5, $6
     }' system.gro > lipid_only.gro

# Produce separate .gro files for each lipid type	 
	 
awk '{	if ( $1 ~ /POPC/ )   
        printf "%9s %5s %5s %7.3f %7.3f %7.3f\n", $1, $2, $3, $4, $5, $6
     }' lipid_only.gro > popc_only.gro

awk '{	if ( $1 ~ /POPE/ )   
        printf "%9s %5s %5s %7.3f %7.3f %7.3f\n", $1, $2, $3, $4, $5, $6
     }' lipid_only.gro > pope_only.gro

awk '{	if ( $1 ~ /POPS/ )   
        printf "%9s %5s %5s %7.3f %7.3f %7.3f\n", $1, $2, $3, $4, $5, $6
     }' lipid_only.gro > pops_only.gro

awk '{	if ( $1 ~ /POPA/ )   
        printf "%9s %5s %5s %7.3f %7.3f %7.3f\n", $1, $2, $3, $4, $5, $6
     }' lipid_only.gro > popa_only.gro

awk '{	if ( $1 ~ /PSM/ )   
        printf "%9s %5s %5s %7.3f %7.3f %7.3f\n", $1, $2, $3, $4, $5, $6
     }' lipid_only.gro > psm_only.gro		

awk '{	if ( $1 ~ /CHL1/ )   
        printf "%9s %5s %5s %7.3f %7.3f %7.3f\n", $1, $2, $3, $4, $5, $6
     }' lipid_only.gro > CHL_only.gro	 
	 

# Start splitting head and tails
	 
echo "### Starting to split heads and tails..."	 

# POPC heads

awk '
BEGIN{
  # variables
  order_count=28
  order_array[1]="N"
  order_array[2]="C12"
  order_array[3]="C13"
  order_array[4]="C14"
  order_array[5]="C15"
  order_array[6]="H12A"
  order_array[7]="H12B"
  order_array[8]="H13A"
  order_array[9]="H13B"
  order_array[10]="H13C"
  order_array[11]="H14A"
  order_array[12]="H14B"
  order_array[13]="H14C"
  order_array[14]="H15A"
  order_array[15]="H15B"
  order_array[16]="H15C"
  order_array[17]="C11"
  order_array[18]="H11A"
  order_array[19]="H11B"
  order_array[20]="P"
  order_array[21]="O13"
  order_array[22]="O14"
  order_array[23]="O12"
  order_array[24]="O11"
  order_array[25]="C1"
  order_array[26]="HA"
  order_array[27]="HB"
  order_array[28]="C2"

  # constants
  count=1
  lipid_count=1
}
{
  particle_array[count] = $0
  particle_type_array[count] = $2
  if( $2 == "N")
  {
    lipid_array[lipid_count] = count
    lipid_count += 1
  }
  count+=1
}
END{
   # count already incremented in line 26
   # lipid_count already incremented in line 24
  lipid_array[lipid_count] = count

  for(i=1; i<lipid_count; i++) # i is index of lipid_array
  {
    for(j=1; j<=order_count; j++) # j is index of order_array
    {
      for(k=lipid_array[i]; k<lipid_array[i+1]; k++) # k in index of particle_array
      {
        if( particle_type_array[k] == order_array[j] )
          print particle_array[k]
      }
    }
  }
}
' popc_only.gro > popc_heads.gro	 
	 
# POPC tails
	 
awk '
BEGIN{
# variables
order_count=108
   order_array[1]="C2"
   order_array[2]="C1" 
   order_array[3]="HS"
   order_array[4]="O21"
   order_array[5]="C21"
   order_array[6]="O22"
   order_array[7]="C22"
   order_array[8]="H2R"
   order_array[9]="H2S"
  order_array[10]="C3"
  order_array[11]="HX"
  order_array[12]="HY"
  order_array[13]="O31"
  order_array[14]="C31"
  order_array[15]="O32"
  order_array[16]="C32"
  order_array[17]="H2X"
  order_array[18]="H2Y"
  order_array[19]="C23"
  order_array[20]="H3R"
  order_array[21]="H3S"
  order_array[22]="C24"
  order_array[23]="H4R"
  order_array[24]="H4S"
  order_array[25]="C25"
  order_array[26]="H5R"
  order_array[27]="H5S"
  order_array[28]="C26"
  order_array[29]="H6R"
  order_array[30]="H6S"
  order_array[31]="C27"
  order_array[32]="H7R"
  order_array[33]="H7S"
  order_array[34]="C28"
  order_array[35]="H8R"
  order_array[36]="H8S"
  order_array[37]="C29"
  order_array[38]="H91"
  order_array[39]="C210"
  order_array[40]="H101"
  order_array[41]="C211"
  order_array[42]="H11R"
  order_array[43]="H11S"
  order_array[44]="C212"
  order_array[45]="H12R"
  order_array[46]="H12S"
  order_array[47]="C213"
  order_array[48]="H13R"
  order_array[49]="H13S"
  order_array[50]="C214"
  order_array[51]="H14R"
  order_array[52]="H14S"
  order_array[53]="C215"
  order_array[54]="H15R"
  order_array[55]="H15S"
  order_array[56]="C216"
  order_array[57]="H16R"
  order_array[58]="H16S"
  order_array[59]="C217"
  order_array[60]="H17R"
  order_array[61]="H17S"
  order_array[62]="C218"
  order_array[63]="H18R"
  order_array[64]="H18S"
  order_array[65]="H18T"
  order_array[66]="C33"
  order_array[67]="H3X"
  order_array[68]="H3Y"
  order_array[69]="C34"
  order_array[70]="H4X"
  order_array[71]="H4Y"
  order_array[72]="C35"
  order_array[73]="H5X"
  order_array[74]="H5Y"
  order_array[75]="C36"
  order_array[76]="H6X"
  order_array[77]="H6Y"
  order_array[78]="C37"
  order_array[79]="H7X"
  order_array[80]="H7Y"
  order_array[81]="C38"
  order_array[82]="H8X"
  order_array[83]="H8Y"
  order_array[84]="C39"
  order_array[85]="H9X"
  order_array[86]="H9Y"
  order_array[87]="C310"
  order_array[88]="H10X"
  order_array[89]="H10Y"
  order_array[90]="C311"
  order_array[91]="H11X"
  order_array[92]="H11Y"
  order_array[93]="C312"
  order_array[94]="H12X"
  order_array[95]="H12Y"
  order_array[96]="C313"
  order_array[97]="H13X"
  order_array[98]="H13Y"
  order_array[99]="C314"
 order_array[100]="H14X"
 order_array[101]="H14Y"
 order_array[102]="C315"
 order_array[103]="H15X"
 order_array[104]="H15Y"
 order_array[105]="C316"
 order_array[106]="H16X"
 order_array[107]="H16Y"
 order_array[108]="H16Z"

  # constants
  count=1
  lipid_count=1
}
{
  particle_array[count] = $0
  particle_type_array[count] = $2
  if( $2 == "N")
  {
    lipid_array[lipid_count] = count
    lipid_count += 1
  }
  count+=1
}
END{


  lipid_array[lipid_count] = count

  for(i=1; i<lipid_count; i++) # i is index of lipid_array
  {
    for(j=1; j<=order_count; j++) # j is index of order_array
    {
      for(k=lipid_array[i]; k<lipid_array[i+1]; k++) # k in index of particle_array
      {
        if( particle_type_array[k] == order_array[j] )
          print particle_array[k]
      }
    }
  }
}
' popc_only.gro > popc_tails.gro

# POPE heads

awk '
BEGIN{
  # variables
  order_count=19
  order_array[1]="N"
  order_array[2]="HN1"
  order_array[3]="HN2"
  order_array[4]="HN3"
  order_array[5]="C12"
  order_array[6]="H12A"
  order_array[7]="H12B"
  order_array[8]="C11"
  order_array[9]="H11A"
  order_array[10]="H11B"
  order_array[11]="P"
  order_array[12]="O13"
  order_array[13]="O14"
  order_array[14]="O11"
  order_array[15]="O12"
  order_array[16]="C1"
  order_array[17]="HA"
  order_array[18]="HB"
  order_array[19]="C2"

  # constants
  count=1
  lipid_count=1
}
{
  particle_array[count] = $0
  particle_type_array[count] = $2
  if( $2 == "N")
  {
    lipid_array[lipid_count] = count
    lipid_count += 1
  }
  count+=1
}
END{
   # count already incremented in line 26
   # lipid_count already incremented in line 24
  lipid_array[lipid_count] = count

  for(i=1; i<lipid_count; i++) # i is index of lipid_array
  {
    for(j=1; j<=order_count; j++) # j is index of order_array
    {
      for(k=lipid_array[i]; k<lipid_array[i+1]; k++) # k in index of particle_array
      {
        if( particle_type_array[k] == order_array[j] )
          print particle_array[k]
      }
    }
  }
}
' pope_only.gro > pope_heads.gro

# POPE tails

awk '
BEGIN{
# variables
order_count=108
   order_array[1]="C2"
   order_array[2]="C1" 
   order_array[3]="HS"
   order_array[4]="O21"
   order_array[5]="C21"
   order_array[6]="O22"
   order_array[7]="C22"
   order_array[8]="H2R"
   order_array[9]="H2S"
  order_array[10]="C3"
  order_array[11]="HX"
  order_array[12]="HY"
  order_array[13]="O31"
  order_array[14]="C31"
  order_array[15]="O32"
  order_array[16]="C32"
  order_array[17]="H2X"
  order_array[18]="H2Y"
  order_array[19]="C23"
  order_array[20]="H3R"
  order_array[21]="H3S"
  order_array[22]="C24"
  order_array[23]="H4R"
  order_array[24]="H4S"
  order_array[25]="C25"
  order_array[26]="H5R"
  order_array[27]="H5S"
  order_array[28]="C26"
  order_array[29]="H6R"
  order_array[30]="H6S"
  order_array[31]="C27"
  order_array[32]="H7R"
  order_array[33]="H7S"
  order_array[34]="C28"
  order_array[35]="H8R"
  order_array[36]="H8S"
  order_array[37]="C29"
  order_array[38]="H91"
  order_array[39]="C210"
  order_array[40]="H101"
  order_array[41]="C211"
  order_array[42]="H11R"
  order_array[43]="H11S"
  order_array[44]="C212"
  order_array[45]="H12R"
  order_array[46]="H12S"
  order_array[47]="C213"
  order_array[48]="H13R"
  order_array[49]="H13S"
  order_array[50]="C214"
  order_array[51]="H14R"
  order_array[52]="H14S"
  order_array[53]="C215"
  order_array[54]="H15R"
  order_array[55]="H15S"
  order_array[56]="C216"
  order_array[57]="H16R"
  order_array[58]="H16S"
  order_array[59]="C217"
  order_array[60]="H17R"
  order_array[61]="H17S"
  order_array[62]="C218"
  order_array[63]="H18R"
  order_array[64]="H18S"
  order_array[65]="H18T"
  order_array[66]="C33"
  order_array[67]="H3X"
  order_array[68]="H3Y"
  order_array[69]="C34"
  order_array[70]="H4X"
  order_array[71]="H4Y"
  order_array[72]="C35"
  order_array[73]="H5X"
  order_array[74]="H5Y"
  order_array[75]="C36"
  order_array[76]="H6X"
  order_array[77]="H6Y"
  order_array[78]="C37"
  order_array[79]="H7X"
  order_array[80]="H7Y"
  order_array[81]="C38"
  order_array[82]="H8X"
  order_array[83]="H8Y"
  order_array[84]="C39"
  order_array[85]="H9X"
  order_array[86]="H9Y"
  order_array[87]="C310"
  order_array[88]="H10X"
  order_array[89]="H10Y"
  order_array[90]="C311"
  order_array[91]="H11X"
  order_array[92]="H11Y"
  order_array[93]="C312"
  order_array[94]="H12X"
  order_array[95]="H12Y"
  order_array[96]="C313"
  order_array[97]="H13X"
  order_array[98]="H13Y"
  order_array[99]="C314"
 order_array[100]="H14X"
 order_array[101]="H14Y"
 order_array[102]="C315"
 order_array[103]="H15X"
 order_array[104]="H15Y"
 order_array[105]="C316"
 order_array[106]="H16X"
 order_array[107]="H16Y"
 order_array[108]="H16Z"

  # constants
  count=1
  lipid_count=1
}
{
  particle_array[count] = $0
  particle_type_array[count] = $2
  if( $2 == "N")
  {
    lipid_array[lipid_count] = count
    lipid_count += 1
  }
  count+=1
}
END{


  lipid_array[lipid_count] = count

  for(i=1; i<lipid_count; i++) # i is index of lipid_array
  {
    for(j=1; j<=order_count; j++) # j is index of order_array
    {
      for(k=lipid_array[i]; k<lipid_array[i+1]; k++) # k in index of particle_array
      {
        if( particle_type_array[k] == order_array[j] )
          print particle_array[k]
      }
    }
  }
}
' pope_only.gro > pope_tails.gro

# POPS heads

awk '
BEGIN{
  # variables
  order_count=21
  order_array[1]="N"
  order_array[2]="HN1"
  order_array[3]="HN2"
  order_array[4]="HN3"
  order_array[5]="C12"
  order_array[6]="H12A"
  order_array[7]="C13"
  order_array[8]="O13A"
  order_array[9]="O13B"
  order_array[10]="C11"
  order_array[11]="H11A"
  order_array[12]="H11B"
  order_array[13]="P"
  order_array[14]="O13"
  order_array[15]="O14"
  order_array[16]="O12"
  order_array[17]="O11"
  order_array[18]="C1"
  order_array[19]="HA"
  order_array[20]="HB"
  order_array[21]="C2"

  # constants
  count=1
  lipid_count=1
}
{
  particle_array[count] = $0
  particle_type_array[count] = $2
  if( $2 == "N")
  {
    lipid_array[lipid_count] = count
    lipid_count += 1
  }
  count+=1
}
END{
   # count already incremented in line 26
   # lipid_count already incremented in line 24
  lipid_array[lipid_count] = count

  for(i=1; i<lipid_count; i++) # i is index of lipid_array
  {
    for(j=1; j<=order_count; j++) # j is index of order_array
    {
      for(k=lipid_array[i]; k<lipid_array[i+1]; k++) # k in index of particle_array
      {
        if( particle_type_array[k] == order_array[j] )
          print particle_array[k]
      }
    }
  }
}
' pops_only.gro > pops_heads.gro

# POPS tails

awk '
BEGIN{
# variables
order_count=108
   order_array[1]="C2"
   order_array[2]="C1" 
   order_array[3]="HS"
   order_array[4]="O21"
   order_array[5]="C21"
   order_array[6]="O22"
   order_array[7]="C22"
   order_array[8]="H2R"
   order_array[9]="H2S"
  order_array[10]="C3"
  order_array[11]="HX"
  order_array[12]="HY"
  order_array[13]="O31"
  order_array[14]="C31"
  order_array[15]="O32"
  order_array[16]="C32"
  order_array[17]="H2X"
  order_array[18]="H2Y"
  order_array[19]="C23"
  order_array[20]="H3R"
  order_array[21]="H3S"
  order_array[22]="C24"
  order_array[23]="H4R"
  order_array[24]="H4S"
  order_array[25]="C25"
  order_array[26]="H5R"
  order_array[27]="H5S"
  order_array[28]="C26"
  order_array[29]="H6R"
  order_array[30]="H6S"
  order_array[31]="C27"
  order_array[32]="H7R"
  order_array[33]="H7S"
  order_array[34]="C28"
  order_array[35]="H8R"
  order_array[36]="H8S"
  order_array[37]="C29"
  order_array[38]="H91"
  order_array[39]="C210"
  order_array[40]="H101"
  order_array[41]="C211"
  order_array[42]="H11R"
  order_array[43]="H11S"
  order_array[44]="C212"
  order_array[45]="H12R"
  order_array[46]="H12S"
  order_array[47]="C213"
  order_array[48]="H13R"
  order_array[49]="H13S"
  order_array[50]="C214"
  order_array[51]="H14R"
  order_array[52]="H14S"
  order_array[53]="C215"
  order_array[54]="H15R"
  order_array[55]="H15S"
  order_array[56]="C216"
  order_array[57]="H16R"
  order_array[58]="H16S"
  order_array[59]="C217"
  order_array[60]="H17R"
  order_array[61]="H17S"
  order_array[62]="C218"
  order_array[63]="H18R"
  order_array[64]="H18S"
  order_array[65]="H18T"
  order_array[66]="C33"
  order_array[67]="H3X"
  order_array[68]="H3Y"
  order_array[69]="C34"
  order_array[70]="H4X"
  order_array[71]="H4Y"
  order_array[72]="C35"
  order_array[73]="H5X"
  order_array[74]="H5Y"
  order_array[75]="C36"
  order_array[76]="H6X"
  order_array[77]="H6Y"
  order_array[78]="C37"
  order_array[79]="H7X"
  order_array[80]="H7Y"
  order_array[81]="C38"
  order_array[82]="H8X"
  order_array[83]="H8Y"
  order_array[84]="C39"
  order_array[85]="H9X"
  order_array[86]="H9Y"
  order_array[87]="C310"
  order_array[88]="H10X"
  order_array[89]="H10Y"
  order_array[90]="C311"
  order_array[91]="H11X"
  order_array[92]="H11Y"
  order_array[93]="C312"
  order_array[94]="H12X"
  order_array[95]="H12Y"
  order_array[96]="C313"
  order_array[97]="H13X"
  order_array[98]="H13Y"
  order_array[99]="C314"
 order_array[100]="H14X"
 order_array[101]="H14Y"
 order_array[102]="C315"
 order_array[103]="H15X"
 order_array[104]="H15Y"
 order_array[105]="C316"
 order_array[106]="H16X"
 order_array[107]="H16Y"
 order_array[108]="H16Z"

  # constants
  count=1
  lipid_count=1
}
{
  particle_array[count] = $0
  particle_type_array[count] = $2
  if( $2 == "N")
  {
    lipid_array[lipid_count] = count
    lipid_count += 1
  }
  count+=1
}
END{


  lipid_array[lipid_count] = count

  for(i=1; i<lipid_count; i++) # i is index of lipid_array
  {
    for(j=1; j<=order_count; j++) # j is index of order_array
    {
      for(k=lipid_array[i]; k<lipid_array[i+1]; k++) # k in index of particle_array
      {
        if( particle_type_array[k] == order_array[j] )
          print particle_array[k]
      }
    }
  }
}
' pops_only.gro > pops_tails.gro

# POPA heads

awk '
BEGIN{
  # variables
  order_count=10
  order_array[1]="P"
  order_array[2]="O13"
  order_array[3]="O14"
  order_array[4]="O12"
  order_array[5]="H12"
  order_array[6]="O11"
  order_array[7]="C1"
  order_array[8]="HA"
  order_array[9]="HB"
  order_array[10]="C2"

  # constants
  count=1
  lipid_count=1
}
{
  particle_array[count] = $0
  particle_type_array[count] = $2
  if( $2 == "P")
  {
    lipid_array[lipid_count] = count
    lipid_count += 1
  }
  count+=1
}
END{
   # count already incremented in line 26
   # lipid_count already incremented in line 24
  lipid_array[lipid_count] = count

  for(i=1; i<lipid_count; i++) # i is index of lipid_array
  {
    for(j=1; j<=order_count; j++) # j is index of order_array
    {
      for(k=lipid_array[i]; k<lipid_array[i+1]; k++) # k in index of particle_array
      {
        if( particle_type_array[k] == order_array[j] )
          print particle_array[k]
      }
    }
  }
}
' popa_only.gro > popa_heads.gro

# POPA tails

awk '
BEGIN{
# variables
order_count=108
  order_array[1]="C2"
  order_array[2]="C1" 
  order_array[3]="HS"
  order_array[4]="O21"
  order_array[5]="C21"
  order_array[6]="O22"
  order_array[7]="C22"
  order_array[8]="H2R"
  order_array[9]="H2S"
  order_array[10]="C3"
  order_array[11]="HX"
  order_array[12]="HY"
  order_array[13]="O31"
  order_array[14]="C31"
  order_array[15]="O32"
  order_array[16]="C32"
  order_array[17]="H2X"
  order_array[18]="H2Y"
  order_array[19]="C23"
  order_array[20]="H3R"
  order_array[21]="H3S"
  order_array[22]="C24"
  order_array[23]="H4R"
  order_array[24]="H4S"
  order_array[25]="C25"
  order_array[26]="H5R"
  order_array[27]="H5S"
  order_array[28]="C26"
  order_array[29]="H6R"
  order_array[30]="H6S"
  order_array[31]="C27"
  order_array[32]="H7R"
  order_array[33]="H7S"
  order_array[34]="C28"
  order_array[35]="H8R"
  order_array[36]="H8S"
  order_array[37]="C29"
  order_array[38]="H91"
  order_array[39]="C210"
  order_array[40]="H101"
  order_array[41]="C211"
  order_array[42]="H11R"
  order_array[43]="H11S"
  order_array[44]="C212"
  order_array[45]="H12R"
  order_array[46]="H12S"
  order_array[47]="C213"
  order_array[48]="H13R"
  order_array[49]="H13S"
  order_array[50]="C214"
  order_array[51]="H14R"
  order_array[52]="H14S"
  order_array[53]="C215"
  order_array[54]="H15R"
  order_array[55]="H15S"
  order_array[56]="C216"
  order_array[57]="H16R"
  order_array[58]="H16S"
  order_array[59]="C217"
  order_array[60]="H17R"
  order_array[61]="H17S"
  order_array[62]="C218"
  order_array[63]="H18R"
  order_array[64]="H18S"
  order_array[65]="H18T"
  order_array[66]="C33"
  order_array[67]="H3X"
  order_array[68]="H3Y"
  order_array[69]="C34"
  order_array[70]="H4X"
  order_array[71]="H4Y"
  order_array[72]="C35"
  order_array[73]="H5X"
  order_array[74]="H5Y"
  order_array[75]="C36"
  order_array[76]="H6X"
  order_array[77]="H6Y"
  order_array[78]="C37"
  order_array[79]="H7X"
  order_array[80]="H7Y"
  order_array[81]="C38"
  order_array[82]="H8X"
  order_array[83]="H8Y"
  order_array[84]="C39"
  order_array[85]="H9X"
  order_array[86]="H9Y"
  order_array[87]="C310"
  order_array[88]="H10X"
  order_array[89]="H10Y"
  order_array[90]="C311"
  order_array[91]="H11X"
  order_array[92]="H11Y"
  order_array[93]="C312"
  order_array[94]="H12X"
  order_array[95]="H12Y"
  order_array[96]="C313"
  order_array[97]="H13X"
  order_array[98]="H13Y"
  order_array[99]="C314"
  order_array[100]="H14X"
  order_array[101]="H14Y"
  order_array[102]="C315"
  order_array[103]="H15X"
  order_array[104]="H15Y"
  order_array[105]="C316"
  order_array[106]="H16X"
  order_array[107]="H16Y"
  order_array[108]="H16Z"

  # constants
  count=1
  lipid_count=1
}
{
  particle_array[count] = $0
  particle_type_array[count] = $2
  if( $2 == "P")
  {
    lipid_array[lipid_count] = count
    lipid_count += 1
  }
  count+=1
}
END{


  lipid_array[lipid_count] = count

  for(i=1; i<lipid_count; i++) # i is index of lipid_array
  {
    for(j=1; j<=order_count; j++) # j is index of order_array
    {
      for(k=lipid_array[i]; k<lipid_array[i+1]; k++) # k in index of particle_array
      {
        if( particle_type_array[k] == order_array[j] )
          print particle_array[k]
      }
    }
  }
}
' popa_only.gro > popa_tails.gro

# PSM renaming

awk '{ if ( $2 == "C1S")
		{ $2=" C1" }
		if ( $2 == "H1S")
			{ $2 = " HA" }
		if ( $2 == "H1T")
			{ $2 = " HB" }
		if ( $2 == "C2S")
			{ $2 = " C2" }
		if ( $2 == "H2S")
			{ $2 = " HS" }
		if ( $2 == "C3S")
			{ $2 = "C21" }
		if ( $2 == "H3S")
			{ $2 = "H21" }
		if ( $2 == "O3")
			{ $2 = "OH1" }
		if ( $2 == "HO3")
			{ $2 = "HO1" }
		if ( $2 == "C4S")
			{ $2 = "C22" }
		if ( $2 == "H4S")
			{ $2 = "H22" }
		if ( $2 == "C5S")
			{ $2 = "C23" }
		if ( $2 == "H5S")
			{ $2 = "H23" }
		if ( $2 == "C6S")
			{ $2 = "C24" }
		if ( $2 == "H6S")
			{ $2 = "H4S" }	
		if ( $2 == "H6T")
			{ $2 = "H4R" }	
		if ( $2 == "C7S")
			{ $2 = "C25" }
		if ( $2 == "H7S")
			{ $2 = "H5S" }	
		if ( $2 == "H7T")
			{ $2 = "H5R" }	
		if ( $2 == "C8S")
			{ $2 = "C26" }
		if ( $2 == "H8S")
			{ $2 = "H6S" }	
		if ( $2 == "H8T")
			{ $2 = "H6R" }
		if ( $2 == "C9S")
			{ $2 = "C27" }
		if ( $2 == "H9S")
			{ $2 = "H7S" }	
		if ( $2 == "H9T")
			{ $2 = "H7R" }
		if ( $2 == "C10S")
			{ $2 = " C28" }
		if ( $2 == "H10S")
			{ $2 = " H8S" }	
		if ( $2 == "H10T")
			{ $2 = " H8R" }
		if ( $2 == "C11S")
			{ $2 = "C29" }
		if ( $2 == "H11S")
			{ $2 = " H9S" }	
		if ( $2 == "H11T")
			{ $2 = " H9R" }		
		if ( $2 == "C12S")
			{ $2 = "C210" }
		if ( $2 == "H12S")
			{ $2 = "H10S" }	
		if ( $2 == "H12T")
			{ $2 = "H10R" }
		if ( $2 == "C13S")
			{ $2 = "C211" }
		if ( $2 == "H13S")
			{ $2 = "H11S" }	
		if ( $2 == "H13T")
			{ $2 = "H11R" }
		if ( $2 == "C14S")
			{ $2 = "C212" }
		if ( $2 == "H14S")
			{ $2 = "H12S" }	
		if ( $2 == "H14T")
			{ $2 = "H12R" }
		if ( $2 == "C15S")
			{ $2 = "C213" }
		if ( $2 == "H15S")
			{ $2 = "H13S" }	
		if ( $2 == "H15T")
			{ $2 = "H13R" }	
		if ( $2 == "C16S")
			{ $2 = "C214" }
		if ( $2 == "H16S")
			{ $2 = "H14S" }	
		if ( $2 == "H16T")
			{ $2 = "H14R" }
		if ( $2 == "C17S")
			{ $2 = "C215" }
		if ( $2 == "H17S")
			{ $2 = "H15S" }	
		if ( $2 == "H17T")
			{ $2 = "H15R" }	
		if ( $2 == "C18S")
			{ $2 = "C216" }
		if ( $2 == "H18S")
			{ $2 = "H16S" }	
		if ( $2 == "H18T")
			{ $2 = "H16R" }		
		if ( $2 == "H18U")
			{ $2 = "H16T" }		
		if ( $2 == "NF")
			{ $2 = "N1" }
		if ( $2 == "HNF")
			{ $2 = " H1" }
		if ( $2 == "C1F")
			{ $2 = "C31" }
		if ( $2 == "OF")
			{ $2 = "O32" }
		if ( $2 == "C2F")
			{ $2 = "C32" }
		if ( $2 == "H2F")
			{ $2 = "H2X" }
		if ( $2 == "H2G")
			{ $2 = "H2Y" }
		if ( $2 == "C3F")
			{ $2 = "C33" }
		if ( $2 == "H3F")
			{ $2 = "H3X" }
		if ( $2 == "H3G")
			{ $2 = "H3Y" }	
		if ( $2 == "C4F")
			{ $2 = "C34" }
		if ( $2 == "H4F")
			{ $2 = "H4X" }
		if ( $2 == "H4G")
			{ $2 = "H4Y" }
		if ( $2 == "C5F")
			{ $2 = "C35" }
		if ( $2 == "H5F")
			{ $2 = "H5X" }
		if ( $2 == "H5G")
			{ $2 = "H5Y" }
		if ( $2 == "C6F")
			{ $2 = "C36" }
		if ( $2 == "H6F")
			{ $2 = "H6X" }
		if ( $2 == "H6G")
			{ $2 = "H6Y" }
		if ( $2 == "C7F")
			{ $2 = "C37" }
		if ( $2 == "H7F")
			{ $2 = "H7X" }
		if ( $2 == "H7G")
			{ $2 = "H7Y" }
		if ( $2 == "C8F")
			{ $2 = "C38" }
		if ( $2 == "H8F")
			{ $2 = "H8X" }
		if ( $2 == "H8G")
			{ $2 = "H8Y" }
		if ( $2 == "C9F")
			{ $2 = "C39" }
		if ( $2 == "H9F")
			{ $2 = "H9X" }
		if ( $2 == "H9G")
			{ $2 = "H9Y" }
		if ( $2 == "C10F")
			{ $2 = "C310" }
		if ( $2 == "H10F")
			{ $2 = "H10X" }
		if ( $2 == "H10G")
			{ $2 = "H10Y" }	
		if ( $2 == "C11F")
			{ $2 = "C311" }
		if ( $2 == "H11F")
			{ $2 = "H11X" }
		if ( $2 == "H11G")
			{ $2 = "H11Y" }
		if ( $2 == "C12F")
			{ $2 = "C312" }
		if ( $2 == "H12F")
			{ $2 = "H12X" }
		if ( $2 == "H12G")
			{ $2 = "H12Y" }
		if ( $2 == "C13F")
			{ $2 = "C313" }
		if ( $2 == "H13F")
			{ $2 = "H13X" }
		if ( $2 == "H13G")
			{ $2 = "H13Y" }	
		if ( $2 == "C14F")
			{ $2 = "C314" }
		if ( $2 == "H14F")
			{ $2 = "H14X" }
		if ( $2 == "H14G")
			{ $2 = "H14Y" }
		if ( $2 == "C15F")
			{ $2 = "C315" }
		if ( $2 == "H15F")
			{ $2 = "H15X" }
		if ( $2 == "H15G")
			{ $2 = "H15Y" }
		if ( $2 == "C16F")
			{ $2 = "C316" }
		if ( $2 == "H16F")
			{ $2 = "H16X" }
		if ( $2 == "H16G")
			{ $2 = "H16Y" }				
		if ( $2 == "H16H")
			{ $2 = "H16Z" }		
           printf "%8s %6s %5s %7.3f %7.3f %7.3f\n", $1, $2, $3, $4, $5, $6
     }' psm_only.gro > psm_newnames.gro


# PSM heads

awk '
BEGIN{
  # variables
  order_count=28
  order_array[1]="N"
  order_array[2]="C13"
  order_array[3]="H13A"
  order_array[4]="H13B"
  order_array[5]="H13C"
  order_array[6]="C14"
  order_array[7]="H14A"
  order_array[8]="H14B"
  order_array[9]="H14C"
  order_array[10]="C15"
  order_array[11]="H15A"
  order_array[12]="H15B"
  order_array[13]="H15C"
  order_array[14]="C12"
  order_array[15]="H12A"
  order_array[16]="H12B"
  order_array[17]="C11"
  order_array[18]="H11A"
  order_array[19]="H11B"
  order_array[20]="P"
  order_array[21]="O13"
  order_array[22]="O14"
  order_array[23]="O11"
  order_array[24]="O12"
  order_array[25]="C1"
  order_array[26]="HA"
  order_array[27]="HB"
  order_array[28]="C2"

  # constants
  count=1
  lipid_count=1
}
{
  particle_array[count] = $0
  particle_type_array[count] = $2
  if( $2 == "N")
  {
    lipid_array[lipid_count] = count
    lipid_count += 1
  }
  count+=1
}
END{
   # count already incremented in line 26
   # lipid_count already incremented in line 24
  lipid_array[lipid_count] = count

  for(i=1; i<lipid_count; i++) # i is index of lipid_array
  {
    for(j=1; j<=order_count; j++) # j is index of order_array
    {
      for(k=lipid_array[i]; k<lipid_array[i+1]; k++) # k in index of particle_array
      {
        if( particle_type_array[k] == order_array[j] )
          print particle_array[k]
      }
    }
  }
}
' psm_newnames.gro > psm_heads.gro

# PSM tails

awk '
BEGIN{
# variables
order_count=101
  order_array[1]="C2"
  order_array[2]="C1"
  order_array[3]="HS"
  order_array[4]="N1"
  order_array[5]="H1"
  order_array[6]="C31"
  order_array[7]="O32"
  order_array[8]="C32"
  order_array[9]="H2X"
  order_array[10]="H2Y"
  order_array[11]="C33"
  order_array[12]="H3X"
  order_array[13]="H3Y"
  order_array[14]="C34"
  order_array[15]="H4X"
  order_array[16]="H4Y"
  order_array[17]="C35"
  order_array[18]="H5X"
  order_array[19]="H5Y"
  order_array[20]="C36"
  order_array[21]="H6X"
  order_array[22]="H6Y"
  order_array[23]="C37"
  order_array[24]="H7X"
  order_array[25]="H7Y"
  order_array[26]="C38"
  order_array[27]="H8X"
  order_array[28]="H8Y"
  order_array[29]="C39"
  order_array[30]="H9X"
  order_array[31]="H9Y"
  order_array[32]="C310"
  order_array[33]="H10X"
  order_array[34]="H10Y"
  order_array[35]="C311"
  order_array[36]="H11X"
  order_array[37]="H11Y"
  order_array[38]="C312"
  order_array[39]="H12X"
  order_array[40]="H12Y"
  order_array[41]="C313"
  order_array[42]="H13X"
  order_array[43]="H13Y"
  order_array[44]="C314"
  order_array[45]="H14X"
  order_array[46]="H14Y"
  order_array[47]="C315"
  order_array[48]="H15X"
  order_array[49]="H15Y"
  order_array[50]="C316"
  order_array[51]="H16X"
  order_array[52]="H16Y"
  order_array[53]="H16Z"
  order_array[54]="C21"
  order_array[55]="H21"
  order_array[56]="OH1"
  order_array[57]="HO1"
  order_array[58]="C22"
  order_array[59]="H22"
  order_array[60]="C23"
  order_array[61]="H23"
  order_array[62]="C24"
  order_array[63]="H4S"
  order_array[64]="H4R"
  order_array[65]="C25"
  order_array[66]="H5S"
  order_array[67]="H5R"
  order_array[68]="C26"
  order_array[69]="H6S"
  order_array[70]="H6R"
  order_array[71]="C27"
  order_array[72]="H7S"
  order_array[73]="H7R"
  order_array[74]="C28"
  order_array[75]="H8S"
  order_array[76]="H8R"
  order_array[77]="C29"
  order_array[78]="H9S"
  order_array[79]="H9R"
  order_array[80]="C210"
  order_array[81]="H10S"
  order_array[82]="H10R"
  order_array[83]="C211"
  order_array[84]="H11S"
  order_array[85]="H11R"
  order_array[86]="C212"
  order_array[87]="H12S"
  order_array[88]="H12R"
  order_array[89]="C213"
  order_array[90]="H13S"
  order_array[91]="H13R"
  order_array[92]="C214"
  order_array[93]="H14S"
  order_array[94]="H14R"
  order_array[95]="C215"
  order_array[96]="H15S"
  order_array[97]="H15R"
  order_array[98]="C216"
  order_array[99]="H16S"
  order_array[100]="H16R"
  order_array[101]="H16T"

  # constants
  count=1
  lipid_count=1
}
{
  particle_array[count] = $0
  particle_type_array[count] = $2
  if( $2 == "N")
  {
    lipid_array[lipid_count] = count
    lipid_count += 1
  }
  count+=1
}
END{


  lipid_array[lipid_count] = count

  for(i=1; i<lipid_count; i++) # i is index of lipid_array
  {
    for(j=1; j<=order_count; j++) # j is index of order_array
    {
      for(k=lipid_array[i]; k<lipid_array[i+1]; k++) # k in index of particle_array
      {
        if( particle_type_array[k] == order_array[j] )
          print particle_array[k]
      }
    }
  }
}
' psm_newnames.gro > psm_tails.gro


# CHOLESTEROL


awk '
BEGIN{
# variables
order_count=74
   order_array[1]="C3"
   order_array[2]="O3"
   order_array[3]="H3'\''"
   order_array[4]="H3"
   order_array[5]="C4"
   order_array[6]="H4A"
   order_array[7]="H4B"
   order_array[8]="C5"
   order_array[9]="C6"
  order_array[10]="H6"
  order_array[11]="C7"
  order_array[12]="H7A"
  order_array[13]="H7B"
  order_array[14]="C8"
  order_array[15]="H8"
  order_array[16]="C14"
  order_array[17]="H14"
  order_array[18]="C15"
  order_array[19]="H15A"
  order_array[20]="H15B"
  order_array[21]="C16"
  order_array[22]="H16A"
  order_array[23]="H16B"
  order_array[24]="C17"
  order_array[25]="H17"
  order_array[26]="C13"
  order_array[27]="C18"
  order_array[28]="H18A"
  order_array[29]="H18B"
  order_array[30]="H18C"
  order_array[31]="C12"
  order_array[32]="H12A"
  order_array[33]="H12B"
  order_array[34]="C11"
  order_array[35]="H11A"
  order_array[36]="H11B"
  order_array[37]="C9"
  order_array[38]="H9"
  order_array[39]="C10"
  order_array[40]="C19"
  order_array[41]="H19A"
  order_array[42]="H19B"
  order_array[43]="H19C"
  order_array[44]="C1"
  order_array[45]="H1A"
  order_array[46]="H1B"
  order_array[47]="C2"
  order_array[48]="H2A"
  order_array[49]="H2B"
  order_array[50]="C20"
  order_array[51]="H20"
  order_array[52]="C21"
  order_array[53]="H21A"
  order_array[54]="H21B"
  order_array[55]="H21C"
  order_array[56]="C22"
  order_array[57]="H22A"
  order_array[58]="H22B"
  order_array[59]="C23"
  order_array[60]="H23A"
  order_array[61]="H23B"
  order_array[62]="C24"
  order_array[63]="H24A"
  order_array[64]="H24B"
  order_array[65]="C25"
  order_array[66]="H25"
  order_array[67]="C26"
  order_array[68]="H26A"
  order_array[69]="H26B"
  order_array[70]="H26C"
  order_array[71]="C27"
  order_array[72]="H27A"
  order_array[73]="H27B"
  order_array[74]="H27C"

  # constants
  count=1
  lipid_count=1
}
{
  particle_array[count] = $0
  particle_type_array[count] = $2
  if( $2 == "C3")
  {
    lipid_array[lipid_count] = count
    lipid_count += 1
  }
  count+=1
}
END{


  lipid_array[lipid_count] = count

  for(i=1; i<lipid_count; i++) # i is index of lipid_array
  {
    for(j=1; j<=order_count; j++) # j is index of order_array
    {
      for(k=lipid_array[i]; k<lipid_array[i+1]; k++) # k in index of particle_array
      {
        if( particle_type_array[k] == order_array[j] )
          print particle_array[k]
      }
    }
  }
}
' CHL_only.gro > cholesterol.gro

# Rename the residues and adding virtual sites

echo "### Renaming the residues and adding virtual sites..."

awk '{sub("POPC","PC  "); print}' popc_heads.gro > file.tmp && mv file.tmp popc_heads.gro
awk '{sub("C2 ","DH "); print}' popc_heads.gro > file.tmp && mv file.tmp popc_heads.gro
awk '{sub("POPC","PO  "); print}' popc_tails.gro > file.tmp && mv file.tmp popc_tails.gro
awk '{sub(" C1 ","DDH "); print}' popc_tails.gro > file.tmp && mv file.tmp popc_tails.gro
awk '{sub("POPE","PE  "); print}' pope_heads.gro > file.tmp && mv file.tmp pope_heads.gro
awk '{sub("C2 ","DH "); print}' pope_heads.gro > file.tmp && mv file.tmp pope_heads.gro
awk '{sub("POPE","PO  "); print}' pope_tails.gro > file.tmp && mv file.tmp pope_tails.gro
awk '{sub(" C1 ","DDH "); print}' pope_tails.gro > file.tmp && mv file.tmp pope_tails.gro
awk '{sub("POPS","PS  "); print}' pops_heads.gro > file.tmp && mv file.tmp pops_heads.gro
awk '{sub("C2 ","DH "); print}' pops_heads.gro > file.tmp && mv file.tmp pops_heads.gro
awk '{sub("POPS","PO  "); print}' pops_tails.gro > file.tmp && mv file.tmp pops_tails.gro
awk '{sub(" C1 ","DDH "); print}' pops_tails.gro > file.tmp && mv file.tmp pops_tails.gro
awk '{sub("POPA","PA  "); print}' popa_heads.gro > file.tmp && mv file.tmp popa_heads.gro
awk '{sub("C2 ","DH "); print}' popa_heads.gro > file.tmp && mv file.tmp popa_heads.gro
awk '{sub("POPA","PO  "); print}' popa_tails.gro > file.tmp && mv file.tmp popa_tails.gro
awk '{sub(" C1 ","DDH "); print}' popa_tails.gro > file.tmp && mv file.tmp popa_tails.gro
awk '{sub("PSM","SM "); print}' psm_heads.gro > file.tmp && mv file.tmp psm_heads.gro
awk '{sub("C2 ","DH "); print}' psm_heads.gro > file.tmp && mv file.tmp psm_heads.gro
awk '{sub("PSM ","PO16"); print}' psm_tails.gro > file.tmp && mv file.tmp psm_tails.gro
awk '{sub(" C1 ","DDH "); print}' psm_tails.gro > file.tmp && mv file.tmp psm_tails.gro

# Merge all files together
	 
cat pope_heads.gro popc_heads.gro psm_heads.gro pops_heads.gro popa_heads.gro > heads.gro
cat pope_tails.gro popc_tails.gro pops_tails.gro popa_tails.gro psm_tails.gro > tails.gro
cat heads.gro tails.gro cholesterol.gro > tmp.gro

# Get the number of atoms and box size

awk 'END{ print NR }' tmp.gro > header.txt
tail -n 1 step7_10.gro > footer.txt

# Fix the .gro file format ; please check which command works well with your file and uncomment it out.

awk 'BEGIN{FIELDWIDTHS="10 5 8 8 8 8"}{sub(/^ /,"",$3);sub(/ $/,"",$2);print $1 $2 $3 $4 $5 $6}' tmp.gro > system_split.gro

#awk '{print substr($0,1,15) substr($0,17)}' tmp.gro > system_split.gro

#awk -v FIELDWIDTHS='15 1 *' '{print $1 $3}' tmp.gro > system_split.gro

#rev tmp.gro | 
#awk '
#  match($0,/^\S+(\s+\S+){3}/){
#    val1=substr($0,RSTART,RLENGTH)
#    val2=substr($0,RSTART+RLENGTH)
#    sub(/^  /,"",val2)
#    $0=val1 val2
#  }
#1
#' | rev > system_split.gro

cat header.txt system_split.gro footer.txt > tmp_file.gro

# Set a title for the .gro file

awk -F\| 'BEGIN {print "Split Membrane"};{print};END {print ""}' tmp_file.gro > tmp_split_lipid_only.gro

# Fix the box size

X=`awk '{print $1}' footer.txt`
Y=`awk '{print $2}' footer.txt`
Z=8.5
gmx editconf -f tmp_split_lipid_only.gro -o split_lipid_only.gro -box $X $Y $Z -c

echo
echo "The split membrane .gro file (split_lipid_only.gro) is ready!"
echo

# Make the topol.top file; you need to have the initial_topol.top file in your current folder

#cp initial_topol.top topol_noW.top

var1=`grep -c "PE" tmp.gro`
var2=19
mol_PE=`echo $((var1 / var2))`

var3=`grep -c "PC" tmp.gro`
var4=28
mol_PC=`echo $((var3 / var4))`

var5=`grep -c "PS" tmp.gro`
var6=21
mol_PS=`echo $((var5 / var6))`

var7=`grep -c "PA" tmp.gro`
var8=10
mol_PA=`echo $((var7 / var8))`

var9=`grep -c "SM" tmp.gro`
mol_SM=`echo $((var9 / var4))`

var11=`grep -c "PO " tmp.gro`
var12=108
mol_PO=`echo $((var11 / var12))`

var13=`grep -c "PO16" tmp.gro`
var14=101
mol_PO16=`echo $((var13 / var14))`

var15=`grep -c "CHL1" tmp.gro`
var16=74
mol_CHL1=`echo $((var15 / var16))`

# Get lipid types

lipids=`tail -n +3 split_lipid_only.gro | head -n-1 | awk -v FIELDWIDTHS="5 4 99" '{print $2}' | uniq` 
lipid_actual=`echo $lipids | awk -v lipid_list=$list '
BEGIN{
  lipid_array[1]="PE"
  lipid_array[2]="PC"
  lipid_array[3]="SM"
  lipid_array[4]="PS"
  lipid_array[5]="PA"
  lipid_array[6]="PO"
  lipid_array[7]="PO16"
  lipid_array[8]="CHL1"
  lipid_count=8
}
{
  for(i=1; i<=NF; i++) 
    for(j=1; j<=lipid_count; j++ )
      if($i == lipid_array[j])
        print $i
}'` 

printf '%s\n' $lipid_actual > lipid_types.top
cat initial_topol.top lipid_types.top > topol_noW.top

for i in $lipid_actual; do
	var=`echo $((mol_$i))`
    gawk -i inplace -v mol="$var" -v lip="$i" '{ if ( $1 == lip )
        { $10=mol }
    print
	}' topol_noW.top
done

# Solvate the system

gmx solvate -cp split_lipid_only.gro -cs spc216.gro -o ref_W.gro -p topol_noW.top

perl water_deletor.pl -in ref_W.gro -out system_sol_fix.gro -ref O32 -middle C316 -nwater 3 -v [yes/no]\n > number_of_waters.txt

NUM_=`grep "water molecules remain" number_of_waters.txt`
NUM=`echo $NUM_ | awk -F" " '{print $1}'`
gawk -i inplace -v NUM="$NUM" '{ if ( $1 == "SOL" )
        { $2=NUM }
    print
	}' topol_noW.top

cp topol_noW.top topol_W.top
mv \#topol_noW.top*  topol_noW.top

# Add ions, please calculate the correct number of ions based on the number of water molecules in the system

touch ions.mdp
gmx grompp -f ions.mdp -c system_sol_fix.gro -p topol_W.top -r system_sol_fix.gro -o ions.tpr |& tee -a info.txt
grep "System has non-zero total charge" info.txt

output=`grep "System has non-zero total charge:" info.txt`
NET_CHARGE=`echo $output | awk '{print $NF}'`

echo "please calculate the correct number of ions based on the number of water molecules in the system {Hint: ($NUM)*0.002772"}
echo -n "Write the Number of CL atoms: "
read N
N="$N"
echo -n "Write the Number of NA atoms {Hint: Number of CL atoms - (0$NET_CHARGE)}: "
read P
P="$P"
echo a SOL | gmx_mpi genion -s ions.tpr -o box.gro -p topol_W.top -pname NA -nname CL -np $P -nn $N

#echo a SOL | gmx_mpi genion -s ions.tpr -o box.gro -p topol_W.top -pname NA -nname CL -neutral -conc 0.154

mv topol_W.top topol.top

echo
echo "#############################"
echo "### Making Reference File ###"
echo "#############################"
echo

# Make reference file (ref.gro)

awk -F '|' -v OFS='|' '{sub(/^.............../, "& "); print}' box.gro > box_.gro
head -n-1 box_.gro > file.tmp && mv file.tmp box_.gro
tail -n+3 box_.gro > file.tmp && mv file.tmp box_.gro

gmx editconf -f step7_10.gro -o tmp_step7_10.gro -box $X $Y $Z -c

echo SYSTEM | gmx traj -f tmp_step7_10.gro -s step7_10.tpr -n index.ndx -ox system.xvg -com
com_system=`awk '{print $NF}' system.xvg | tail -n 1` 

awk -v var="$com_system" '{ if ( $2 == "OW" )
        { $6=var }
       printf "%9s %5s %5s %7.3f %7.3f %7.3f\n", $1, $2, $3, $4, $5, $6
	 }' box_.gro > tmp_ref.gro  

awk -v var="$com_system" '{ if ( $1 ~ /PE/ && $2 == "C1" )
        { $6=var }
       printf "%9s %5s %5s %7.3f %7.3f %7.3f\n", $1, $2, $3, $4, $5, $6
	 }' tmp_ref.gro > file.tmp && mv file.tmp tmp_ref.gro

awk -v var="$com_system" '{ if ( $1 ~ /PC/ && $2 == "C1" )
        { $6=var }
       printf "%9s %5s %5s %7.3f %7.3f %7.3f\n", $1, $2, $3, $4, $5, $6
	 }' tmp_ref.gro > file.tmp && mv file.tmp tmp_ref.gro

awk -v var="$com_system" '{ if ( $1 ~ /PS/ && $2 == "C1" )
        { $6=var }
       printf "%9s %5s %5s %7.3f %7.3f %7.3f\n", $1, $2, $3, $4, $5, $6
	 }' tmp_ref.gro > file.tmp && mv file.tmp tmp_ref.gro	 

awk -v var="$com_system" '{ if ( $1 ~ /PA/ && $2 == "C1" )
        { $6=var }
       printf "%9s %5s %5s %7.3f %7.3f %7.3f\n", $1, $2, $3, $4, $5, $6
	 }' tmp_ref.gro > file.tmp && mv file.tmp tmp_ref.gro

awk -v var="$com_system" '{ if ( $1 ~ /SM/ && $2 == "C1" )
        { $6=var }
       printf "%9s %5s %5s %7.3f %7.3f %7.3f\n", $1, $2, $3, $4, $5, $6
	 }' tmp_ref.gro > file.tmp && mv file.tmp tmp_ref.gro	 
	 	 
#awk -v var="$com_system" '{ if ( $1 ~ /PO/ && $2 == "C316" )
#        { $6=var }
#      printf "%9s %5s %5s %7.3f %7.3f %7.3f\n", $1, $2, $3, $4, $5, $6
#	 }' tmp_ref.gro > file.tmp && mv file.tmp tmp_ref.gro
#
#awk -v var="$com_system" '{ if ( $1 ~ /PO/ && $2 == "C218" )
#        { $6=var }
#      printf "%9s %5s %5s %7.3f %7.3f %7.3f\n", $1, $2, $3, $4, $5, $6
#	 }' tmp_ref.gro > file.tmp && mv file.tmp tmp_ref.gro
#
#awk -v var="$com_system" '{ if ( $1 ~ /PO16/ && $2 == "C216" )
#        { $6=var }
#      printf "%9s %5s %5s %7.3f %7.3f %7.3f\n", $1, $2, $3, $4, $5, $6
#	 }' tmp_ref.gro > file.tmp && mv file.tmp tmp_ref.gro		 	 	 
	 
# Get the number of atoms and box size
	 
awk 'END{ print NR }' box_.gro > header.txt
tail -n 1 box.gro > footer.txt

# Fix the .gro file format

awk 'BEGIN{FIELDWIDTHS="5 5 5 8 8 8 8"}{sub(/^ /,"",$4);sub(/ $/,"",$3);print $1 $2 $3 $4 $5 $6 $7}' tmp_ref.gro > ref_.gro	 
	 
cat header.txt ref_.gro footer.txt > ref.gro

# Set a title for the .gro file

awk -F\| 'BEGIN {print "Reference File"};{print};END {print ""}' ref.gro > file.tmp && mv file.tmp ref.gro

# Remove extra files

rm p*.gro
rm system_split.gro
rm heads.gro
rm tails.gro
rm system.gro
rm lipid*
rm system_sol_fix.gro
rm *.txt
rm tmp*
rm mdout*
rm *.xvg
rm box_.gro
rm ref_.gro
rm ref_W.gro
rm topol_noW.top
rm ions.*
rm CHL*
rm chol*
rm \#*

# Make index file

residues=`tail -n +3 box.gro | head -n-1 | awk -v FIELDWIDTHS="5 4 99" '{print $2}' | uniq` 
residue_actual=`echo $residues | awk -v lipid_list=$list '
BEGIN{
  lipid_array[1]="PE"
  lipid_array[2]="PC"
  lipid_array[3]="SM"
  lipid_array[4]="PS"
  lipid_array[5]="PA"
  lipid_array[6]="PO"
  lipid_array[7]="PO16"
  lipid_array[8]="CHL1"
  lipid_count=8
}
{
  for(i=1; i<=NF; i++) 
    for(j=1; j<=lipid_count; j++ )
      if($i == lipid_array[j])
        print $i
}'`

    {
echo "keep 0"
for i in $residue_actual
do
  echo -n "r $i |"
done
echo
echo "name 1 MEMB"
echo "r SOL | r NA | r CL"
echo "name 2 SOLV"
echo q
    } | gmx make_ndx -f box.gro -o groups.ndx

# Run a short minimization

gmx grompp -f EM.mdp -c box.gro -n groups.ndx -r ref.gro -p topol.top -o split_bilayer.tpr
gmx_mpi mdrun -v -deffnm split_bilayer

echo
echo "###########################################"
echo "###  The split bilayer is ready to use. ###"
echo "###        Enjoy simulation!            ###"
echo "###########################################"
echo
