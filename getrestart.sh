#!/bin/tcsh -f
#
# get some restart files from HPSS or specified directory
# using the specified basename, find the file of the form basename????.????.u
# with the newest time stamp (????.????)
#
#  getrestart.sh basname HPSS [uvw,h5,vor]        search HPSS
#  getrestart.sh basname path [uvw,h5,vor]        search path 
#  getrestart.sh basname path [uvw,h5,vor]  all    search both
#
#  $1 = basename of restart file. 
#  $2 = HPSS   get file off of HPSS
#     = path   get file from specified path
#  $3   uvw or h5 or vor.  default = uvw
#  $4   "all"
#
# status=0  success
# status=1  failed
#

set PSIGET = "psi get"
#set PSIGET = "echo psi get disabled: "

#set PSIVM = "psi mv"
set PSIMV  = "echo psi mv disabled: "

set name = $1
set fpath = $2
set ext = w

if ($#argv >= 3 ) then
   if ($3 == uvw ) then
      set ext = w
   else if ($3 == h5) then
      set ext = h5
   else if ($3 == vor) then
      set ext = vor
   endif
endif

set searchall = no
if ($#argv >= 4 ) then
   set searchall = $4
endif


if ($searchall == all) then
   #search $fpath for newest restart file
   echo searching both HPSS and $fpath     
   set resnamels = `\ls {$fpath}/{$name}*.{$ext} | sort | tail -1`
   if ($resnamels =="") then
      echo {$fpath}/{$name}
      echo "no restart files in $fpath"
      set resnamels = aaaaaaaaaaa
   endif
   #search PSI
   set resnamepsi = `psi ls dns/{$name}/{$name}\*.{$ext} | sort | tail -1`
   if ($resnamepsi =="") then
      echo "no restart files on HPSS"
      set resnamepsi = aaaaaaaaaaa
   endif
   set resnamels = `basename $resnamels`
   set resnamepsi = `basename $resnamepsi`
   rm -f tempsort
   echo $resnamels > tempsort
   echo $resnamepsi >> tempsort 
   set resnamew=`sort tempsort | tail -1`
   rm -f tempsort
   if ($resnamew == $resnamels ) then
      echo "using restart files in directory"
      #fpath already set
   else if ($resnamew == $resnamepsi ) then
      echo "using restart files from HPSS"
      set fpath = HPSS
   else
      echo "no restart files found"
      exit 1
   endif
endif


if ($fpath == HPSS) then

   #search HPSS for newest restart file
   set resnamew = `psi ls dns/{$name}/{$name}\*.{$ext} | sort | tail -1`
   if ($resnamew =="") then
      echo "Error finding restart file.  Exit"
      exit 1
   else
      echo "Using HPSS restart files: " 
      set nametime = `basename $resnamew .{$ext}`
      if ( $ext == "w" ) then
         #set resnameu = `psi ls  dns/{$name}/{$nametime}\*.u | sort | tail -1`
         #set resnamev = `psi ls  dns/{$name}/{$nametime}\*.v | sort | tail -1`
         set resnameu = dns/{$name}/{$nametime}.u
         set resnamev = dns/{$name}/{$nametime}.v
         echo $resnameu
         echo $resnamev
      endif
      if ( $ext == "vor" ) then
         set resnamet = dns/{$name}/{$nametime}.tracer 
         echo $resnamet
      endif
      echo $resnamew
   endif
   if ( $ext == "w" ) then
      set resnameu2 = `basename $resnameu`
      set resnamev2 = `basename $resnamev`
   endif
   if ( $ext == "vor" ) then
      set resnamet2 = `basename $resnamet`
   endif
   set resnamew2 = `basename $resnamew`

   # check to see if files are left over from last run: 
   \rm -f restart.*

   if ( $ext == "w" ) then
      if !(-e $resnameu2) then
         $PSIGET $resnameu
      endif
      if !(-e $resnamev2) then
         $PSIGET $resnamev
      endif
      if !(-e $resnamew2) then
         $PSIGET $resnamew
      endif
      \ln -s $resnameu2  restart.u
      \ln -s $resnamev2  restart.v
      \ln -s $resnamew2  restart.{$ext}
      if !(-e restart.u) then
         echo "No restart.u file"
         $PSIMV $resnameu $resnameu.bak
         $PSIMV $resnamev $resnamev.bak
         $PSIMV $resnamew $resnamew.bak
         exit 1
      endif
   else  if ( $ext == "h5" ) then
      if !(-e $resnamew2) then
         $PSIGET $resnamew
      endif
      \ln -s $resnamew2  restart.{$ext}
      if !(-e restart.{$ext}) then
         echo "No restart.{$ext} file"
         $PSIMV $resnamew $resnamew.bak
         exit 1
      endif
   else
      if !(-e $resnamet2) then
         $PSIGET $resnamet
      endif
      if !(-e $resnamew2) then
         $PSIGET $resnamew
      endif
      \ln -s $resnamet2  restart.tracer
      \ln -s $resnamew2  restart.{$ext}
      if !(-e restart.tracer) then
         echo "No restart.tracer file"
         $PSIMV $resnamet $resnamet.bak
         $PSIMV $resnamew $resnamew.bak
         exit 1
      endif
   endif

else

   #search $fpath for newest restart file
   set resnamew = `ls {$fpath}/{$name}*.{$ext} | sort | tail -1`
   if ($resnamew =="") then
      echo {$fpath}/{$name}
      echo "Error finding restart file.  Exit"
      exit 1
   else
      echo "Using directory restart files: " 
      set nametime = `basename $resnamew .{$ext}`
      if ( $ext == "w" ) then
         #set resnameu = `ls  {$fpath}/{$nametime}*.u | sort | tail -1`
         #set resnamev = `ls  {$fpath}/{$nametime}*.v | sort | tail -1`
         set resnameu = {$fpath}/{$nametime}.u 
         set resnamev = {$fpath}/{$nametime}.v 
         echo $resnameu
         echo $resnamev
      endif
      if ( $ext == "vor" ) then
         set resnamet = {$fpath}/{$nametime}.tracer 
         echo $resnamet
      endif
      echo $resnamew
   endif

   \rm -f restart.*
   if ( $ext == "w" ) then
      \ln -s $resnameu  restart.u
      \ln -s $resnamev  restart.v
      \ln -s $resnamew  restart.{$ext}
      if !(-e restart.{$ext}) then
         echo "No restart.{$ext} file"
         # move them out of the way in case this time is corrupt:
         # then next run will pick up earlier backup:
         mv $resnameu $resnameu.bak
         mv $resnamev $resnamev.bak
         mv $resnamew $resnamew.bak
         exit 1
      endif
   else if ( $ext == "h5" ) then
      \ln -s $resnamew  restart.{$ext}
      if !(-e restart.{$ext}) then
         echo "No restart.{$ext} file"
         # move them out of the way in case this time is corrupt:
         # then next run will pick up earlier backup:
         mv $resnamew $resnamew.bak
         exit 1
       endif
   else
      \ln -s $resnamet  restart.tracer
      \ln -s $resnamew  restart.{$ext}
      if !(-e restart.{$ext}) then
         echo "No restart.{$ext} file"
         # move them out of the way in case this time is corrupt:
         # then next run will pick up earlier backup:
         mv $resnamet $resnamet.bak
         mv $resnamew $resnamew.bak
         exit 1
      endif
   endif
endif

exit 0











