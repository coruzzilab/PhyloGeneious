# *** Perl ***
# OrthologID phylogenetic tree class
#
# Copyright (C) 2006-2011 Ernest K. Lee
#
#    This file is part of OrthologID.
#
#    OrthologID is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    OrthologID is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with OrthologID.  If not, see <http://www.gnu.org/licenses/>.
#
# Author: Ernest K Lee <elee@amnh.org>
#
package OrthologID::PhyloTree;
use strict;
our ( @ISA, @EXPORT );

use Exporter;
@ISA    = qw( Exporter );
@EXPORT = qw( new printNodes orthologGroups );

my $unknownSp = "Unknown";

#
# Constructor
# Args: parenthetical tree (string), ortholog species list (ref)
#
sub new($$) {
	my $package       = shift;
	my $parenTree     = shift;
	my $speciesPrefix = shift;
	my $treeObj       = toTreeObj( $parenTree, $speciesPrefix );
	bless $treeObj, $package;
	return $treeObj;
}

#
# Generate a tree object from a parenthetical tree
#
sub toTreeObj {
	my $pTree    = shift;
	my $spPrefix = shift;

	my $nodeNum = 0;

	# Create root node; remove root parentheses
	my $treeObj = newNode( undef, $nodeNum++ );
	my $currNode;

	#my $currNode = $treeObj;
	#$pTree =~ s/^\((.*)\)$/$1/;

	# Traverse the parenthetical tree
	while ( $pTree =~ /([\(\),])([^\(\),]*)/g ) {
		my ( $paren, $name ) = ( $1, $2 );
		if ( $paren eq "(" ) {
			if ( !defined($currNode) ) {
				$currNode = $treeObj;
			}
			else {
				my $newNode = newNode( $currNode, $nodeNum++ );
				push( @{ $currNode->{"child"} }, $newNode );
				$currNode = $newNode;
			}
		}
		elsif ( $paren eq ")" ) {

			# Calculate evo status of current node
			my %spMem;
			my $spDup = 0; # number of species appear in more than 1 children
			foreach my $ch ( @{ $currNode->{"child"} } ) {
				$currNode->{"eStatus"} *= $ch->{"eStatus"};
				last if $currNode->{"eStatus"} == 0;  # nothing to check
				foreach (@$spPrefix) {
					if ( defined( $ch->{"spMember"}->{$_} ) ) {
						if ( defined( $spMem{$_} ) ) {

							# Same species in more than one child
							#$currNode->{"eStatus"} = 0;
							#last;
							$spDup++;
						}
						else {
							$spMem{$_} = 1;
						}
					}
				}
			}
			# if more than 1 species present in all children, treat as duplication
			$currNode->{"eStatus"} = 0 if $spDup > 1;
				

			# Pure duplications if only one species
			$currNode->{"eStatus"} = 1
			  if keys( %{ $currNode->{"spMember"} } ) == 1;

			my $parent = $currNode->{"parent"};
			# Done if we are at root
			last if !defined $parent;

			# Propagate species members to parent
			foreach (@$spPrefix) {
				if ( defined( $currNode->{"spMember"}->{$_} ) ) {
					$parent->{"spMember"}->{$_} +=
					  $currNode->{"spMember"}->{$_};
				}
			}

			# Move up
			$currNode = $parent;
		}

		if ( $name ne "" ) {

			# Create leaf child (taxon) and add species member
			my $newNode = newNode( $currNode, $nodeNum++, $name );
			push( @{ $currNode->{"child"} }, $newNode );
			my $spMember = $currNode->{"spMember"};
			my $matched  = 0;
			foreach (@$spPrefix) {
				if ( $name =~ /^$_#/ ) {
					$spMember->{$_}++;    # Add taxon species to count
					$newNode->{"spMember"}->{$_} = 1;
					$matched = 1;
					last;
				}
			}
			if ( !$matched ) {            # should not happen

				# unknown species
				warn "Unknown species ($name) encountered in tree!";
				$spMember->{$unknownSp}++;
				$newNode->{"spMember"}->{$unknownSp}++;
			}
		}
	}
	return $treeObj;
}

#
# Create new tree node.  Returns pointer.
#
sub newNode {
	my ( $parent, $nodeNum, $label ) = @_;
	my %node = (
		num      => $nodeNum,    # Node number in DFS order
		label    => $label,      # Taxon name
		child    => [],
		parent   => $parent,
		spMember => {},          # Taxa species of this node (clade)
		eStatus => 1    # Evolution status of node (pure dup/sp[1] or mix[0])
	);
	return \%node;

}

#
# Print tree nodes and their member species (for debugging)
#
sub printNodes {
	my $treeObj = shift;

	print "node #" . $treeObj->{"num"} . " [" . $treeObj->{"eStatus"} . "]";
	while ( my ( $sp, $num ) = each( %{ $treeObj->{"spMember"} } ) ) {
		print " $sp: $num ";
	}
	print $treeObj->{"label"} if defined $treeObj->{"label"};
	print "\n";
	foreach ( @{ $treeObj->{"child"} } ) {
		printNodes($_);
	}
}


#
# Return list of ortholog groups as list of list refs
# Arguments are list of species of interests
#
sub orthologGroups {
	my $treeNode = shift;
	my @orthSp   = @_;
	my @oGroups  = ();

	if ( $treeNode->{"eStatus"} ) {
		my $spCount = 0;
		foreach (@orthSp) {
			$spCount++ if defined( $treeNode->{"spMember"}->{$_} );
		}
		my @group = ();
		@group = getLeaves( $treeNode, @orthSp ) if $spCount > 1;
		@oGroups = ( \@group ) if @group > 0;
	}
	else {
		foreach my $child ( @{ $treeNode->{"child"} } ) {
			my @groups = orthologGroups( $child, @orthSp );
			push( @oGroups, @groups ) if @groups > 0;
		}
	}
	return @oGroups;

}

#
# Return list of leaves of a node belonging to species in arguments
# getLeaves(@)
#
sub getLeaves {
	my $node   = shift;
	my @sp     = @_;
	my @leaves = ();

	foreach my $child ( @{ $node->{"child"} } ) {
		if ( defined $child->{"label"} ) {
			my $label = $child->{"label"};
			foreach my $sp (@sp) {
				if ( $label =~ /^$sp/ ) {
					push( @leaves, $label );
					last;
				}
			}
		}
		else {
			push( @leaves, getLeaves( $child, @sp ) );
		}
	}

	return @leaves;
}

1;
