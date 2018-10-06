#######################################
# this bash script renames all pdbs   #
# in number format                    #
#######################################

a=1
for i in *.pdb; do
  new=$(printf "%d.pdb" "$a")
  mv -i -- "$i" "$new"
  let a=a+1
done
