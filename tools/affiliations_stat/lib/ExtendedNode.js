/*
 * ExtendedNode 1.0.0 - JavaScript ExtendedNode Library
 *
 * Copyright (c) 2015 INRA
 * 
 * @author: Escudié Frédéric
 * @license: MIT (http://www.opensource.org/licenses/mit-license.php) license.
 */
var ExtendedNode = function( pName, pParent_node, pChild_nodes, pMetadata ) {
    this.name = pName ;
    this.parent = pParent_node ;
	this.children = new Array();
    if( pChild_nodes != null ){
		for( var child_idx = 0 ; child_idx < pChild_nodes.length ; child_idx++ ){
			this.addChild( pChild_nodes[child_idx] );
		}
	}
    this.metadata = new Array();
    if( pMetadata != null ){
		this.metadata = pMetadata
    }
}

ExtendedNode.prototype.hasChild = function( pName ){
	if( pName == null ){
		return this.children.length > 0 ;
	} else {
		var find = false ;
		for( var idx = 0 ; idx < this.children.length && !find ; idx++ ){
			if( pName == this.children[idx] ){
				find = true ;
			}
		}
		return find ;
		
	}
};

ExtendedNode.prototype.getChild = function( pName ){
	if( !this.hasChild( pName ) ){
		throw this.name + " doesn't have child named '" + pName + "'." ;
	}
	return this.children[pName] ;
};

ExtendedNode.prototype.getParent = function(){
	return this.parent ;
};

ExtendedNode.prototype.getNodeByDepth = function( pDepth ){
	var nodes = new Array();
	
	if( pDepth == 0 ){
		nodes.push( this );
	} else if( this.hasChild() ){
		for( var child_idx = 0 ; child_idx < this.children.length ; child_idx++ ){
			var descendant_with_depth = this.children[child_idx].getNodeByDepth( pDepth-1 );
			for( var desc_idx = 0 ; desc_idx < descendant_with_depth.length ; desc_idx++ ){
				nodes.push( descendant_with_depth[desc_idx] );
			}
		}
	}
	
	return nodes ;
};

ExtendedNode.prototype.addChild = function( pChild ){
	if( pChild.parent == null ){
		pChild.parent = this ;
	}
	this.children.push( pChild );
};
		  
ExtendedNode.fromNewick = function( pNewick ){
	var stack = new Array();
	for( var str_idx = 0 ; str_idx < pNewick.length ; str_idx++ ){
		var current_char = pNewick[str_idx] ;
		if( current_char == "(" ){ // Start ancestor node
			stack.push( "(" );
		} else if( current_char == ")" ) { // End ancestor node
			var parent = new ExtendedNode( null, null, null, null );
			while( stack[stack.length-1] != "(" ){
				parent.addChild( stack.pop() );
			}
			stack.pop(); // remove '('
			stack.push( parent );
		} else if( current_char == "," ) {
			// nothing
		} else if( current_char == ";" ) {
			// nothing
		} else if( current_char == " " ) {
			// nothing
		} else if( current_char == ":" ) { // Nodes distance
			if( pNewick[str_idx+1] == "{" ){ // Extended newick (use "{metadata_A: 10, ...}" instead of distance)
				metadata_json = "{" ;
				str_idx++ ;
				nb_open = 1 ;
				nb_closed = 0 ;
				while( nb_open != nb_closed ){
					str_idx++ ;
					metadata_json += pNewick[str_idx] ;
					if( pNewick[str_idx] == "{" ){
						nb_open++ ;
					} else if( pNewick[str_idx] == "}" ){
						nb_closed++ ;
					}
				}
				stack[stack.length-1].metadata = JSON.parse(metadata_json);
			} else { // Standard newick
				var distance = "" ;
				while( pNewick[str_idx+1] != "," && pNewick[str_idx+1] != ")" && pNewick[str_idx+1] != ";" ){
					str_idx++ ;
					distance += pNewick[str_idx] ;
				}
				stack[stack.length-1].metadata['dist'] = distance ;
			}
		} else { // Node name
			var previous_char = null ;
			if( str_idx != 0 ){
				previous_char = pNewick[str_idx -1] ;
			}
			if( current_char == '"' ){
				stop_markers = ['"'];
				str_idx++ ;
			} else {
				stop_markers = [',', ')', ':', ';'];
			}
			var node_name = pNewick[str_idx] ;
			while( stop_markers.indexOf(pNewick[str_idx+1]) == -1 ){
				str_idx++ ;
				node_name += pNewick[str_idx] ;
			}
			if( stop_markers.indexOf('"') != -1 ){
				str_idx++ ;
			}
			if( previous_char == ")" ){
				stack[stack.length-1].name = node_name ;
			} else {
				var node = new ExtendedNode( node_name, null, null, null );
				stack.push( node );
			}
		}
	}
	return stack[0] ;
};
		  
ExtendedNode.prototype.keepOnlySamples = function( kept_samples ){
	var to_remove = false ;
	if( this.hasChild() ){
		var children_to_remove = new Array();
		for( var child_idx = 0 ; child_idx < this.children.length ; child_idx++ ){
			if( this.children[child_idx].keepOnlySamples(kept_samples) ){
				children_to_remove.push( child_idx );
			}
		}
		var offset = 0 ;
		for( var remove_idx = 0 ; remove_idx < children_to_remove.length ; remove_idx++ ){
			this.children.splice( children_to_remove[remove_idx] - offset, 1 );
			offset++ ;
		}
		if( this.children.length == 0 ){
			to_remove = true ;
		}
	} else {
		var keep_metadata = false ;
		var new_metadata = {} ;
		for( var sample_idx = 0 ; sample_idx < kept_samples.length ; sample_idx++ ){
			if( this.metadata.hasOwnProperty(kept_samples[sample_idx]) ){
				new_metadata[kept_samples[sample_idx]] = this.metadata[kept_samples[sample_idx]] ;
				keep_metadata = true ;
			}
		}
		if( keep_metadata ){
			this.metadata = new_metadata ;
		} else {
			to_remove = true ;
		};
	}
	return to_remove ;
};

ExtendedNode.prototype.toJson = function(){
	var node_json = {
		'name': (this.name != null)? this.name : "",
		'metadata': this.metadata
	};
	if( this.hasChild() ){
		var children_json = new Array();
		for( var child_idx = 0 ; child_idx < this.children.length ; child_idx++ ){
			children_json.push( this.children[child_idx].toJson() );
		}
		node_json['children'] = children_json ;
	}
	// Add size
	var metadata_keys = Object.keys( this.metadata );
	if( metadata_keys.length != 0 ){
		var size = 0 ;
		for( var metadata_idx = 0 ; metadata_idx < metadata_keys.length ; metadata_idx++ ){
			size += this.metadata[metadata_keys[metadata_idx]] ;
		}
		node_json['size'] = size ;
	}
	return node_json ;
};