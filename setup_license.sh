#!/bin/bash
NOTNT="You have not agreed to all the terms of the license, or
  password file contains the wrong password

  For info on non-public versions of TNT, contact
  Pablo Goloboff (pablogolo@csnat.unt.edu.ar)

  For further info, check the License page at:

  http://www.lillo.org.ar/phylogeny/tnt/files/LicenseAgreement.htm
  
  OrthologID cannot be used without TNT. Exiting with error."

OID_HOME=`pwd`
sed -i "s|OID_HOME|${OID_HOME}|" setup_rundir.sh
sed -i "s|\(OID_HOME=\).*|\1${OID_HOME}|" testdata/run.sh

echo "Starting PhyloGeneious setup..."
if [ ! -f $OID_HOME/.passwordfile.tnt ]; then
	echo "You are about to go to the TNT license agreement. Please follow the instructions on the screen and enter 'zzz' when you are finished."
	echo "Press enter on the keyboard when you are ready to begin."
	read -p EMPTY
	while [ ! -f $OID_HOME/.passwordfile.tnt ]; do
#		/tnt-linux/TNT-bin/tnt 2>/terminalErr
#		cat /terminalErr | grep -q "Error opening terminal"
#		if [ $? == 0 ]; then 
#			echo "Terminal failed: initiating manual agreement"
			echo "TNT Personal Use License (PUL)

 The TNT PUL allows you to use TNT free of charge for personal
 use only.

 If you do not belong to this category, you will have to purchase
 a commercial license. Do not hesitate to contact the authors of
 the program in this matter.

(press 'y' to accept, 'n' to decline)"
			i=0
			while [ $i -lt 1 ]; do
				read RESPONSE
				if [ "${RESPONSE//[[:space:]]/}" == n ]; then
					echo $NOTNT
					exit 1
				elif [ "${RESPONSE//[[:space:]]/}" == y ]; then
					i=1
				fi
			done
			echo "TNT PUL terms and conditions

 Version 1.6, April 2023

 The Authors grant you the right to use the software product as
 defined in # 1 according to the following provisions. If you do not
 agree to all conditions set forth by this license, you may not use
 the product, because only The Authors as the product's owners can
 give you permission to use it.

 The PUL allows you to download the TNT binaries for personal use,
 but it does not give you the right to redistribute these binaries.
 So you may not put them onto your own websites or other mirrors.

(press 'y' to accept, 'n' to decline)"
			while [ $i -lt 2 ]; do
				read RESPONSE
				if [ "${RESPONSE//[[:space:]]/}" == n ]; then
					echo $NOTNT
					exit 1
				elif [ "${RESPONSE//[[:space:]]/}" == y ]; then
					i=2
				fi
			done
			echo "# 1 Subject of license. "Product", as referred to in this License,
 shall be the binary software package "TNT", which allows for
 calculating and evaluating phylogenetic trees. The Product consists
 of executable files in machine code for the Windows 2000/XP, Linux
 and Mac operating systems as well as other data files as required by
 the executable files at run-time and documentation in electronic
 form. "The Authors", as referred to in this License, shall be Pablo
 A. Goloboff (who can be contacted at pablogolo@csnat.unt.edu.ar),
 James S. Farris (who can be contacted at msl-farr@nrm.se), and
 Kevin C. Nixon (who can be contacted at kcn2@cornell.edu

(press 'y' to accept, 'n' to decline)"
			while [ $i -lt 3 ]; do
				read RESPONSE
				if [ "${RESPONSE//[[:space:]]/}" == n ]; then
					echo $NOTNT
					exit 1
				elif [ "${RESPONSE//[[:space:]]/}" == y ]; then
					i=3
				fi
			done
			echo "# 2 Grant of license. The Authors grant you a personal right to
 install and execute the Product on a Host Computer for Personal Use
 or Educational Use. "Personal Use" requires that you use the product
 on the same Host Computer where you installed it yourself and that no
 more than one client connect to that Host Computer at a time for the
 purpose of running a phylogenetic analysis.

(press 'y' to accept, 'n' to decline)"
			while [ $i -lt 4 ]; do
				read RESPONSE
				if [ "${RESPONSE//[[:space:]]/}" == n ]; then
					echo $NOTNT
					exit 1
				elif [ "${RESPONSE//[[:space:]]/}" == y ]; then
					i=4
				fi
			done
			echo "# 3 Reservation of rights. Any use beyond the provisions of # 1 is
 prohibited. The Authors reserve all copyrights and other intellectual
 property rights. This includes, but is not limited to, the right to
 modify, make available or public, rent out, lease, lend or otherwise
 distribute the Product. This does not apply as far as applicable law
 may require it or The Authors grant you additional rights of use in
 a separate license in writing.

(press 'y' to accept, 'n' to decline)"
			while [ $i -lt 5 ]; do
				read RESPONSE
				if [ "${RESPONSE//[[:space:]]/}" == n ]; then
					echo $NOTNT
					exit 1
				elif [ "${RESPONSE//[[:space:]]/}" == y ]; then
					i=5
				fi
			done
			echo "# 4 Termination. This License shall be valid indefinitely. The Authors
 may terminate the License only for material causes. In particular,
 such a material cause can be a violation of the usage terms or a
 breach of other essential duties from this contract. After termina_
 tion, you are required to delete and destroy all remaining copies of
 the Product. This includes, but is not limited to, installed copies
 and backups.

(press 'y' to accept, 'n' to decline)"
			while [ $i -lt 6 ]; do
				read RESPONSE
				if [ "${RESPONSE//[[:space:]]/}" == n ]; then
					echo $NOTNT
					exit 1
				elif [ "${RESPONSE//[[:space:]]/}" == y ]; then
					i=6
				fi
			done
			echo "# 5 No warranties. Since you have not paid for the use of the Product,
 there is no warranty for it, to the extent permitted by applicable law.
 The Authors provide the Product "as is" without warranty of any kind,
 either expressed or implied, including, but not limited to, the implied
 warranties of merchantability and fitness for a particular purpose.
 The entire risk as to the quality and performance of the Product is
 with you. Should it prove defective, you assume the cost of all
 necessary servicing, repair, or correction.

(press 'y' to accept, 'n' to decline)"
			while [ $i -lt 7 ]; do
				read RESPONSE
				if [ "${RESPONSE//[[:space:]]/}" == n ]; then
					echo $NOTNT
					exit 1
				elif [ "${RESPONSE//[[:space:]]/}" == y ]; then
					i=7
				fi
			done
			echo "# 6 Publication of results.  You are free to publish in a scientific
 journal any results derived from use of the program.  In such a case,
 you have to (1) acknowledge that the program is being made available
 with the sponsorship of the Willi Hennig Society, and (2) cite the
 paper describing the program:

     Goloboff, P., & Morales, M. 2023. TNT version 1.6, with a graphical
       interface for MacOs and Linux, including new routines in parallel.
       Cladistics DOI 10.1111/cla.12524



(press 'y' to accept, 'n' to decline)"
			while [ $i -lt 8 ]; do
				read RESPONSE
				if [ "${RESPONSE//[[:space:]]/}" == n ]; then
					echo $NOTNT
					exit 1
				elif [ "${RESPONSE//[[:space:]]/}" == y ]; then
					i=8
				fi
			done
			echo "# 7 Miscellaneous. There are no license terms beyond the written ones
 in this agreement. Amendments of, additions to and the joint
 revocation of this agreement shall require the written form. The same
 shall apply to the preceding written form requirement. Standard
 business conditions of the parties shall not apply. Place of
 performance and legal venue shall be San Miguel de Tucuman, the
 domicile of Pablo A. Goloboff.  Solely Argentinian law shall apply
 to this agreement.

(press 'y' to accept, 'n' to decline)"
			while [ $i -lt 9 ]; do
				read RESPONSE
				if [ "${RESPONSE//[[:space:]]/}" == n ]; then
					echo $NOTNT
					exit 1
				elif [ "${RESPONSE//[[:space:]]/}" == y ]; then
					i=9
				fi
			done
			echo "Type 'I agree' (yes: all six letters!) to finish...

Do you agree to all the terms and conditions?"
			while [ $i -lt 10 ]; do
				read RESPONSE
				RESPONSE=${RESPONSE^^}
				if [ "${RESPONSE//[[:space:]]/}" == N ]; then
					echo $NOTNT
					exit 1
				elif [ "${RESPONSE//[[:space:]]/}" == "IAGREE" ]; then
					i=10
					#cp $OID_HOME/.passwordfile.tnt /root
					echo palosgomias > $OID_HOME/.passwordfile.tnt
				else
					echo "Type 'I agree' to finish (or 'n' to decline)"
				fi
			done
#		fi
		if [ ! -f $OID_HOME/.passwordfile.tnt ]; then
			echo "TNT license agreement failed. Would you like to try again?"
			echo "(enter 'y' to re-try, or press anything to continue)"
			read RESPONSE
			#echo $RESPONSE
			#echo [ ${RESPONSE} = "Y" ]
			if [ "${RESPONSE//[[:space:]]/}" != [Yy] ]; then
			echo "OrthologID cannot be used without TNT. Exiting container build with error."
			exit 1
			else echo "retry"
			fi
			#./tnt-linux/TNT-bin/tnt
		fi
	done
	echo "Success!"
fi
#echo "Starting pipeline..."
#startoid
