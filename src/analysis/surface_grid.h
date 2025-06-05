#pragma once
#include "spade.h"

#define TOL 8.0
#define PARALLEL_TOL 1e-3

namespace pptb::analysis
{
	enum node_type
	{
		regular_node, // Dual vertex (worth 1 dual face)
		edge_node, // Dual vertex on feature line (worth 2 dual faces)
		vertex_node // Dual vertex on boundary intersection (worth at least 3 dual faces)
	};
	
	struct node_map_t
	{
		node_type type=regular_node;
		int ndual=1;
		std::vector<std::vector<std::size_t>> list; // List of faces forming each dual face (only for edge_node and vertex node types)
	};
	
	template <std::floating_point rtype>
	static auto compute_node_normal(const int idx, const spade::geom::vtk_geom_t<3, 3, rtype>& geom)
	{
		spade::ctrs::array<rtype, 3> nvec = 0.0;

		for (int i = 0; i<geom.node2face_count[idx]; ++i)
		{
			const auto& f = geom.node2face_data[geom.node2face_start[idx]+i];
			nvec += geom.normals[f];
		}

		return nvec / rtype(geom.node2face_count[idx]);
	}

	template <typename geom_t>
	static spade::ctrs::array<std::size_t, 2> find_shared_nodes(const std::size_t f1, const std::size_t f2, const geom_t& geom)
	{
		int node_match = 0;
		spade::ctrs::array<std::size_t, 2> edge_nodes;

		const auto& list1 = geom.faces[f1];
		const auto& list2 = geom.faces[f2];
		for (int i = 0; i<list1.size(); ++i)
		{
			for (int j = 0; j<list2.size(); ++j)
			{
				if (list1[i] == list2[j])
				{
					// Store node
					edge_nodes[node_match] = list1[i];
					++node_match;
				}
			}
		}

		if (node_match != 2)
		{
			print("Shared edge node identification failed!");
			std::cin.get();
		}
		
		return edge_nodes;
	}

	template <typename geom_t>
	static bool edge_is_feature_line(const std::size_t f1, const std::size_t f2, const geom_t& geom)
	{
		using real_t = typename geom_t::value_type;
		using vec_t  = spade::ctrs::array<real_t, 3>;
		
		// My normal
		const vec_t pvec = geom.normals[f1];
		
		// Neighbor normal
		const vec_t qvec = geom.normals[f2];
		
		// Angle between vectors (check for feature line)
		const real_t dot_prod  = spade::utils::max(spade::utils::min(spade::ctrs::dot_prod(pvec, qvec), 1.0), -1.0);
		const real_t angle_loc = acos(dot_prod) * 180.0 / spade::consts::pi;

		return (angle_loc > TOL);
	}
	
	template <typename geom_t>
	static bool shared_node_is_boundary_intersection(const spade::ctrs::array<std::size_t, 2>& list1, const spade::ctrs::array<std::size_t, 2>& list2, const geom_t& geom)
	{
		using real_t = typename geom_t::value_type;
		using vec_t  = spade::ctrs::array<real_t, 3>;
		
		// Vector following edge 1
		vec_t pvec;
		#pragma unroll
		for (int j = 0; j<vec_t::size(); ++j) pvec[j] = geom.points[list1[0]][j] - geom.points[list1[1]][j];
		pvec /= spade::ctrs::array_norm(pvec);

		// Vector following edge 2
		vec_t qvec;
		#pragma unroll
		for (int j = 0; j<vec_t::size(); ++j) qvec[j] = geom.points[list2[1]][j] - geom.points[list2[0]][j];
		qvec /= spade::ctrs::array_norm(qvec);

		// Angle between edges
		const real_t dot_prod  = spade::utils::max(spade::utils::min(spade::ctrs::dot_prod(pvec, qvec), 1.0), -1.0);
		const real_t angle_loc = acos(dot_prod) * 180.0 / spade::consts::pi;

		return (angle_loc > TOL);
	}

	template <typename voro_t, typename geom_t>
	static void reorder_voronoi_nodes(const std::size_t f, voro_t& voronoi, const geom_t& geom)
	{
		using real_t = typename geom_t::value_type;
		using vec_t  = spade::ctrs::array<real_t, 3>;
		using v2_t   = spade::ctrs::array<real_t, 2>;
		using pnt_t  = spade::coords::point_t<real_t>;
		
		// Compute element centroid
		vec_t cen{0.0, 0.0, 0.0};
		for (int p = 0; p<voronoi.face2node_count[f]; ++p)
		{
			const auto& inode = voronoi.face2node_data[voronoi.face2node_start[f]+p];
			cen += voronoi.points[inode];
		}
		cen /= real_t(voronoi.face2node_count[f]);

		// Coordinate system around centroid
		spade::ctrs::array<vec_t, 3> basis;
		
		// Vector from centroid to node 1
		vec_t xyz0;
		spade::ctrs::copy_array(voronoi.points[voronoi.face2node_data[voronoi.face2node_start[f]+0]], xyz0);
		basis[1]  = xyz0 - cen;
		basis[1] /= spade::ctrs::array_norm(basis[1]);

		// Second vector to node 2
		for (int p = 1; p<voronoi.face2node_count[f]; ++p)
		{
			// Get node coordinate
			vec_t xyz;
			spade::ctrs::copy_array(voronoi.points[voronoi.face2node_data[voronoi.face2node_start[f]+p]], xyz);
			vec_t pvec = xyz - cen;
			pvec /= spade::ctrs::array_norm(pvec);

			// Check if these two vectors are aligned
			if ((real_t(1.0) - spade::utils::abs(spade::ctrs::dot_prod(pvec, basis[1]))) > PARALLEL_TOL)
			{
				// Centroid normal (approximate)
				basis[0] = spade::ctrs::cross_prod(basis[1], pvec);
				basis[0] /= spade::ctrs::array_norm(basis[0]);
				break;
			}
		}
		
		// Final basis vector
		basis[2] = spade::ctrs::cross_prod(basis[1], basis[0]);
		basis[2] /= spade::ctrs::array_norm(basis[2]);

		// Sort nodes by computing local in-plane angles w.r.t. centroid
		struct theta_t
		{
			std::size_t id;
			real_t val;
		};
		std::vector<theta_t> theta; theta.resize(voronoi.face2node_count[f], theta_t{0, 0.0});
		for (int p = 0; p<voronoi.face2node_count[f]; ++p)
		{
			// Get node coordinate
			vec_t xyz;
			spade::ctrs::copy_array(voronoi.points[voronoi.face2node_data[voronoi.face2node_start[f]+p]], xyz);
			
			// Vector from centroid to node
			spade::ctrs::array<real_t, 3> pvec = xyz - cen;
			pvec /= spade::ctrs::array_norm(pvec);

			// Distance to plane
			const real_t dist = spade::ctrs::dot_prod(pvec, basis[0]);

			// Project node to plane
			vec_t proj = xyz - dist * basis[0];

			// Vector from centroid to projected node
			proj -= cen;

			// Compute 2D coordinates for point in plane
			v2_t xyz2D;
			xyz2D[0] = spade::ctrs::dot_prod(basis[1], proj);
			xyz2D[1] = spade::ctrs::dot_prod(basis[2], proj);

			// Compute theta
			theta[p].id = voronoi.face2node_data[voronoi.face2node_start[f]+p];
			theta[p].val= atan2(xyz2D[0], xyz2D[1]);
		}

		// Sort thetas
		auto comparator = [&](const auto& a, const auto& b) { return a.val < b.val; };
		std::sort(theta.begin(), theta.end(), comparator);

		// Correct node ordering
		for (int p = 0; p<voronoi.face2node_count[f]; ++p)
			voronoi.face2node_data[voronoi.face2node_start[f]+p] = theta[p].id;
	}

	template <typename voro_t, typename geom_t>
	static void correct_voronoi_normal(const std::size_t f, voro_t& voronoi, const std::size_t fref, const geom_t& geom)
	{
		using real_t = typename geom_t::value_type;
		using vec_t  = spade::ctrs::array<real_t, 3>;
		using pnt_t  = spade::coords::point_t<real_t>;

		// Reference normal
		const auto& ref = geom.normals[fref];

		// Compute element centroid
		vec_t cen{0.0, 0.0, 0.0};
		for (int p = 0; p<voronoi.face2node_count[f]; ++p)
		{
			const auto& inode = voronoi.face2node_data[voronoi.face2node_start[f]+p];
			cen += voronoi.points[inode];
		}
		cen /= real_t(voronoi.face2node_count[f]);
		
		// Vector from centroid to node 1
		vec_t xyz0;
		spade::ctrs::copy_array(voronoi.points[voronoi.face2node_data[voronoi.face2node_start[f]+0]], xyz0);
		vec_t qvec  = xyz0 - cen;
		qvec /= spade::ctrs::array_norm(qvec);

		// Vector from centroid to node 2
		vec_t xyz;
		spade::ctrs::copy_array(voronoi.points[voronoi.face2node_data[voronoi.face2node_start[f]+1]], xyz);
		vec_t pvec = xyz - cen;
		pvec /= spade::ctrs::array_norm(pvec);

		// Cell normal
		vec_t nvec = spade::ctrs::cross_prod(qvec, pvec);
		nvec /= spade::ctrs::array_norm(nvec);

		// Check normal orientation
		const real_t dot_prod = spade::ctrs::dot_prod(nvec, ref);

		// Is this cell inverted?
		if (dot_prod < 0.0)
		{
			std::vector<std::size_t> list; list.resize(voronoi.face2node_count[f], 0);

			// Save list
			#pragma unroll
			for (int p = 0; p<voronoi.face2node_count[f]; ++p)
				list[p] = voronoi.face2node_data[voronoi.face2node_start[f]+p];

			// Invert list
			#pragma unroll
			for (int p = 0; p<voronoi.face2node_count[f]; ++p)
				voronoi.face2node_data[voronoi.face2node_start[f]+p] = list[voronoi.face2node_count[f]-p-1];
		}
	}

	template <typename geom_t, typename map_t, typename voro_t>
	static void restricted_lloyds_algorithm(const int nsmooth, const geom_t& geom, const map_t& map, voro_t& voronoi)
	{
		using real_t = typename geom_t::value_type;
		using vec_t  = spade::ctrs::array<real_t, 3>;
		using pnt_t  = spade::coords::point_t<real_t>;
		
		print("Running restricted lloyds algorithm...");

		// Store seed points
		std::vector<pnt_t> seeds; seeds.resize(geom.points.size());
		for (std::size_t n = 0; n<geom.points.size(); ++n) seeds[n] = geom.points[n];

		// Number of iterations
		const int max_iter = nsmooth;
		int iter           = 0;
		bool complete      = false;
		const real_t RES   = 1e-6;
		real_t error       = 100.0;
		real_t max_err     = 0.0;

		auto compute_residual = [&]()
		{
			// Loop over primal vertices
			max_err            = 0.0;
			error              = 0.0;
			std::size_t npoly  = 0;
			std::size_t nfaces = 0;
			for (std::size_t n = 0; n<map.size(); ++n)
			{
				// Loop over dual faces from this primal vertex
				for (int idual = 0; idual<map[n].ndual; ++idual)
				{
					// Only apply smoothing to regular vertices
					if (map[n].type == regular_node)
					{
						// Count vertices
						++nfaces;

						// Compute dual face centroid
						const auto& dual_centroid = voronoi.centroid(npoly);

						// Get seed location (vertex on primal grid)
						const auto& seed_point    = seeds[n];

						// Compute local error
						real_t err_loc = 0.0;
						#pragma unroll
						for (int i = 0; i<dual_centroid.size(); ++i) err_loc += (dual_centroid[i] - seed_point[i]) * (dual_centroid[i] - seed_point[i]);

						// Accumulate error to global L2 error norm
						error += err_loc;

						// Maximum error
						max_err = (sqrt(err_loc) > max_err) ? sqrt(err_loc) : max_err;
					}

					// Dual face counter
					++npoly;
				}
			}

			// Finalize L2 error norm
			error = sqrt(error / real_t(nfaces));
		};
			
		// Iterative smoother
		while (!complete && iter++<max_iter)
		{			
			// Update seed points
			std::size_t npoly  = 0;
			for (std::size_t n = 0; n<map.size(); ++n)
			{
				// Loop over dual faces from this primal vertex
				for (int idual = 0; idual<map[n].ndual; ++idual)
				{
					// Only apply smoothing to regular vertices
					if (map[n].type == regular_node)
					{
						// Compute dual face centroid
						const auto& dual_centroid = voronoi.centroid(npoly);

						// Update seed location
						seeds[n] = dual_centroid;

						// Maintain surface shape. Project onto underlying geometry
						std::size_t fid;
						seeds[n] = geom.slide_to_closest_point(geom.points[n], seeds[n], 45.0, fid);
					}

					// Dual face counter
					++npoly;
				}
			}

			// Push updated seed points to voronoi structure (update regular dual faces only)
			for (std::size_t f = 0; f<geom.faces.size(); ++f)
			{
				// Compute new centroid on primal mesh
				pnt_t cen{0.0, 0.0, 0.0};
				for (int n = 0; n<geom.num_edges(); ++n)
				{
					const auto& inode = geom.faces[f][n];
					cen[0] += seeds[inode][0];
					cen[1] += seeds[inode][1];
					cen[2] += seeds[inode][2];
				}
				cen /= real_t(geom.num_edges());
				
				// Update Voronoi vertex
				voronoi.points[f] = cen;
			}
			
			// Compute error
			compute_residual();
			print("iter, error, max  = ", iter, error, max_err);
			
			// Termination condition
			//if (error < RES) complete = true;
			if (max_err < RES) complete = true;
		}
	}
	
	template <std::floating_point rtype>
	static auto dual_conversion(const int nsmooth, const spade::geom::vtk_geom_t<3, 3, rtype>& geom)
	{
		// Following --> "Polyhedral Mesh Generation for CFD-Analysis of Complex Structures" by Balafas, Georgios (2014)
		print("Converting to surface grid dual...");
		using real_t = rtype;
		using vec_t  = spade::ctrs::array<real_t, 3>;
		using pnt_t  = spade::coords::point_t<real_t>;

		// Voronoi geometry
		spade::geom::voronoi_geom_t<3, rtype> voronoi;

		// Invalid identifier
		constexpr std::size_t invalid = std::numeric_limits<std::size_t>::max();

		// Step 1: All primal face centroids are marked as dual vertices
		print("   Marking primal centroids...");
		voronoi.points.resize(geom.faces.size());
		for (std::size_t f = 0; f<geom.faces.size(); ++f)
		{
			voronoi.points[f] = geom.centroid(f);
		}

		// Step 2: All primal edges on feature lines are significant. Dual vertices are placed at edge midpoints
		print("   Adding midpoints on primal edges...");
		std::vector<std::size_t> edge_list;
		std::vector<spade::ctrs::array<std::size_t, 2>> pair_list;
		for (std::size_t f = 0; f<geom.faces.size(); ++f)
		{
			// Point centroid
			const auto& xyz = geom.centroid(f);
			
			// Neighbor list
			const auto& neigh_list = geom.edge_neighbors[f];

			// Loop neighbors
			for (int p = 0; p<neigh_list.size(); ++p)
			{
				// Neighbor
				const auto& fneigh = neigh_list[p];

				// Beyond tolerance? Store this edge
				if (edge_is_feature_line(f, fneigh, geom)) edge_list.push_back(spade::utils::pairing_function(f, fneigh));
			}
		}

		// Remove duplicate edges
		std::sort(edge_list.begin(), edge_list.end());
		edge_list.erase(std::unique(edge_list.begin(), edge_list.end()), edge_list.end());

		auto generate_pair_list = [&]()
		{
			// Clear memory
			pair_list.clear();
			
			// Save pairs
			for (std::size_t i = 0; i<edge_list.size(); ++i)
			{
				// Edge ID
				const auto& edgeID = edge_list[i];

				// Faces
				std::size_t f1, f2;
				spade::utils::pairing_function_inv(edgeID, f1, f2);

				// Shared nodes
				auto edge_nodes = find_shared_nodes(f1, f2, geom);

				// Smaller node ID first
				if (edge_nodes[0] > edge_nodes[1])
				{
					std::size_t tmp = edge_nodes[0];
					edge_nodes[0] = edge_nodes[1];
					edge_nodes[1] = tmp;
				}
				pair_list.push_back(edge_nodes);
			}
		};
		generate_pair_list();
		
		// Error checking. Delete any edges where a node is only connected to 1 edge. Erroneous feature line identification
		auto bad_edge = [&](const auto& e)
		{
			// Look up edge
			const auto ipt = spade::algs::binary_search(edge_list, e, invalid);

			// Get node pair
			const auto& edge_nodes = pair_list[ipt];

			// Loop both nodes
			for (int i = 0; i<edge_nodes.size(); ++i)
			{
				// Node ID
				const auto& inode = edge_nodes[i];

				// Loop node face neighbors
				std::vector<std::size_t> connected_edges;
				for (int j = 0; j<geom.node2face_count[inode]; ++j)
				{
					// Face 1
					const auto& face1 = geom.node2face_data[geom.node2face_start[inode]+j];
					
					// Edge neighbors of this face
					for (int k = 0; k<geom.edge_neighbors[face1].size(); ++k)
					{
						// Face 2
						const auto& face2 = geom.edge_neighbors[face1][k];

						// Compute edgeID
						const auto& edgeID = spade::utils::pairing_function(face1, face2);

						// Make sure this is the same edge and that its a feature edge
						const auto ipt2 = spade::algs::binary_search(edge_list, edgeID, invalid);
						
						// Neighbor must be a face node neighbor of our node
						bool good_neighbor = false;
						for (int t = 0; t<geom.node2face_count[inode]; ++t)
							if (face2 == geom.node2face_data[geom.node2face_start[inode]+t])
							{
								good_neighbor = true;
								break;
							}

						// Make sure neighbor is a node face neighbor of originating node
						if (good_neighbor)
						{
							if (ipt2 != invalid)
							{
								connected_edges.push_back(edgeID);
							}
						} // End good neighbor check
					} // End edge neighbor loop
				} // End node face neighbor loop

				// Sort mapped edges and remove duplicates
				std::sort(connected_edges.begin(), connected_edges.end());
				connected_edges.erase(std::unique(connected_edges.begin(), connected_edges.end()), connected_edges.end());
				
				// If we counted no other feature edges, remove this edge
				if (connected_edges.size() == 1) return true;
				
			} // End shared node loop
			
			return false;
		};
		
		bool proceed = true;
		std::size_t curr_size = edge_list.size();
		int max_iter = 1000;
		int iter = 0;
		while (proceed && iter++<max_iter)
		{
			// Find edges to remove
			const auto itt = std::remove_if(edge_list.begin(), edge_list.end(), bad_edge);

			// Delete edges
			edge_list.erase(itt, edge_list.end());

			// Regenerate pair list
			generate_pair_list();

			// New number of edges
			std::size_t new_size = edge_list.size();

			// Exit condition
			if (curr_size == new_size) proceed = false;
			else curr_size = new_size;
		}
		
		// Add edge midpoints to voronoi graph
		for (std::size_t i = 0; i<edge_list.size(); ++i)
		{
			// Shared nodes
			const auto& edge_nodes = pair_list[i];

			// Midpoint coordinate
			const pnt_t midpoint{
				real_t(0.5) * (geom.points[edge_nodes[0]][0] + geom.points[edge_nodes[1]][0]),
				real_t(0.5) * (geom.points[edge_nodes[0]][1] + geom.points[edge_nodes[1]][1]),
				real_t(0.5) * (geom.points[edge_nodes[0]][2] + geom.points[edge_nodes[1]][2])
			};

			// Add to graph
			voronoi.points.push_back(midpoint);
		}

		// Step 3: Primal vertices on feature edge intersections are significant. Add as voronoi vertex
		print("   Marking primal vertices on feature edge intersections...");
		auto idx_comparator = [&](const auto& a, const auto& b)
		{
			if (a[0] != b[0]) return a[0] < b[0];
			if (a[1] != b[1]) return a[1] < b[1];
			return false;
		}; // Sort based on first node, then second
		std::sort(pair_list.begin(), pair_list.end(), idx_comparator);
		std::vector<std::size_t> bndy_int;
		for (std::size_t i = 0; i<pair_list.size()-1; ++i)
		{
			// My edge nodes
			const auto& list1 = pair_list[i];

			// Next edge nodes
			const auto& list2 = pair_list[i+1];

			// Make sure these edges are connected
			if (list1[0] == list2[0])
			{
				// Beyond tolerance? Store this edge
				if (shared_node_is_boundary_intersection(list1, list2, geom)) bndy_int.push_back(list1[0]);
			}
		}

		// Remove duplicates
		std::sort(bndy_int.begin(), bndy_int.end());
		bndy_int.erase(std::unique(bndy_int.begin(), bndy_int.end()), bndy_int.end());
		for (std::size_t i = 0; i<bndy_int.size(); ++i) voronoi.points.push_back(geom.points[bndy_int[i]]);

		// Mapping from primal nodes to dual face type
		std::vector<node_map_t> node_map;
		node_map.resize(geom.points.size());

		// Step 8/9: Form dual faces from dual vertices (steps 4-7 are just marking edges which we will do here too)
		print("   Characterizaing primal vertices...");
		// Set default settings
		for (std::size_t n = 0; n<geom.points.size(); ++n)
		{
			// Get mapping
			auto& map = node_map[n];
			
			// Default --> assume this is a regular vertex
			map.type      = regular_node;
			map.ndual     = 1;
		}

		// Loop feature edges and classify nodes on these edges as edge_node or vertex_node
		std::vector<int> mask; mask.resize(geom.points.size(), 0);
		for (std::size_t e = 0; e<edge_list.size(); ++e)
		{
			// Get edge ID
			const auto& edgeID = edge_list[e];

			// Face connections
			std::size_t f1, f2;
			spade::utils::pairing_function_inv(edgeID, f1, f2);

			// Shared nodes
			const auto& list = find_shared_nodes(f1, f2, geom);
			//print("e, edgeID, f1, f2, edge_nodes = ", e, edgeID, f1, f2, list);

			// Loop shared nodes
			for (int n = 0; n<list.size(); ++n)
			{
				// Node ID
				const auto& inode = list[n];
				//print("   n, inode, xyz              = ", n, inode, geom.points[inode]);
				
				// Have we characterized this node yet?
				std::vector<int> my_pairs;
				if (!mask[inode])
				{
					// Map for this node
					auto& map = node_map[inode];

					// Find how many feature edges we are associated with
					for (std::size_t p = 0; p<pair_list.size(); ++p)
						if (pair_list[p][0] == inode || pair_list[p][1] == inode)
							my_pairs.push_back(p);
					
					// Remove duplicates
					my_pairs.erase(std::unique(my_pairs.begin(), my_pairs.end()), my_pairs.end());
					//print("   my_pairs.size()            = ", my_pairs.size());

					// Set primal vertex type
					if (my_pairs.size() == 2) map.type = edge_node;
					else map.type = vertex_node;

					// Number of dual faces formed by this primal vertex is equal to number of feature edges we're associated with
					map.ndual = my_pairs.size();
					map.list.resize(my_pairs.size());
					std::vector<int> mask_face; mask_face.resize(geom.node2face_count[inode], 0);
					int ipair = 0;
					for (int m = 0; m<geom.node2face_count[inode]; ++m)
					{
						// Reference element is the face we're on
						const auto& fref = geom.node2face_data[geom.node2face_start[inode]+m];
						
						// Get reference normal
						const auto& nvec = geom.normals[fref];
						//print("      m, fref, nvec           = ", m, fref, nvec);
						
						// Have we mapped out this face already?
						if (!mask_face[m]) 
						{
							// Add reference face to new pair list
							map.list[ipair].push_back(fref);

							// Mask this face
							mask_face[m] = 1;
							
							// Loop node face neighbors again and find faces which are on same plane
							for (int k = 0; k<geom.node2face_count[inode]; ++k)
							{
								// Neighbor ID
								const auto& fneigh = geom.node2face_data[geom.node2face_start[inode]+k];

								// Make sure this isn't ourselves and we aren't already mapped
								if (fneigh != fref && !mask_face[k])
								{
									// Neighbor normal
									const auto& pvec = geom.normals[fneigh];

									// Angle between faces
									const real_t dot_prod  = spade::utils::max(spade::utils::min(spade::ctrs::dot_prod(nvec, pvec), 1.0), -1.0);
									const real_t angle_loc = acos(dot_prod) * 180.0 / spade::consts::pi;
									
									//print("         k, fneigh, pvec      = ", k, fneigh, pvec);
									//print("         dot_prod, angle_loc  = ", dot_prod, angle_loc, TOL);

									// Same surface?
									if (angle_loc < TOL)
									{
										// Add to list
										map.list[ipair].push_back(fneigh);
										
										// Mask face
										mask_face[k] = 1;
									}
								}
							}

							// Done with this reference element. Incremenet pair index for prep for next dual face setup
							++ipair;
						}
					}
					
					// Mask this node
					mask[inode] = 1;
				}
			}
		}

		// Allocate memory for storage
		print("   Allocating memory for storage...");
		std::size_t npoly = 0;
		std::size_t nelem = 0;
		for (std::size_t n = 0; n<geom.points.size(); ++n)
		{
			npoly += node_map[n].ndual;
			if (node_map[n].type == regular_node) nelem += geom.node2face_count[n];
			else
			{
				// Loop over dual faces coming from this vertex
				for (int p = 0; p<node_map[n].list.size(); ++p)
				{
					nelem += node_map[n].list[p].size() + 2; // Add two to account for midpoints to close dual face
					if (node_map[n].type == vertex_node) ++nelem; // If were a vertex node, we need to add ourselves to dual vertex count too
				}
			}
		}
		voronoi.face2node_count.resize(npoly);
		voronoi.face2node_start.resize(npoly, 0);
		voronoi.face2node_data.resize(nelem, invalid);
		voronoi.component.resize(npoly, -1);
		
		// Compute and store dual face connectivity data
		print("   Compute and store dual face connectivity...");
		npoly = 0;
		for (std::size_t n = 0; n<geom.points.size(); ++n)
		{
			// Local primal vertex map
			const auto& map = node_map[n];
			
			// Is this a regular node?
			if (map.type == regular_node)
			{ // Form dual face from primal node2face connectivity

				// Number of vertices on dual face
				voronoi.face2node_count[npoly] = geom.node2face_count[n];

				// Start position in vector
				if (npoly > 0) voronoi.face2node_start[npoly] = voronoi.face2node_start[npoly-1] + voronoi.face2node_count[npoly-1];

				// Connectivity storage (node numbers equivalent to primal face ID's since this was counted first)
				for (int p = 0; p<voronoi.face2node_count[npoly]; ++p)
					voronoi.face2node_data[voronoi.face2node_start[npoly]+p] = geom.node2face_data[geom.node2face_start[n]+p];

				// Get component number from any face
				voronoi.component[npoly] = geom.component[geom.node2face_data[geom.node2face_start[n]+0]];

				// Re-order voronoi node ordering. Needs to be sequential and match nvec orientation of original surface
				reorder_voronoi_nodes(npoly, voronoi, geom);

				// Correct voronoi normal
				correct_voronoi_normal(npoly, voronoi, geom.node2face_data[geom.node2face_start[n]+0], geom);
				
				// Count poly faces
				++npoly;
			}
			else
			{ // Form dual face from combination of primal centroids and edge midpoints
					
				// Loop dual faces
				for (int idual = 0; idual<map.ndual; ++idual)
				{
					// Number of vertices on this dual face
					voronoi.face2node_count[npoly] = map.list[idual].size() + 2;
					if (map.type == vertex_node) ++voronoi.face2node_count[npoly];
							
					// Start position in vector
					if (npoly > 0) voronoi.face2node_start[npoly] = voronoi.face2node_start[npoly-1] + voronoi.face2node_count[npoly-1];

					// Get component number from any face
					voronoi.component[npoly] = geom.component[map.list[idual][0]];

					// Nodes connecting this face (don't worry about order yet)
					int nnode = 0;
					for (int p = 0; p<map.list[idual].size(); ++p) // Loop primal faces on associated to this dual face
					{
						// Primal face ID
						const auto& pf = map.list[idual][p];

						// Add primal face to vertex list
						voronoi.face2node_data[voronoi.face2node_start[npoly]+nnode] = pf;
						++nnode;

						// Check neighbors. If a neighbor isn't in list[idual] then it's associated with a different dual face. Mark midpoint and add to list
						for (int k = 0; k<geom.edge_neighbors[pf].size(); ++k)
						{
							// Neighbor ID
							const auto& fneigh = geom.edge_neighbors[pf][k];

							// Neighbor must be a face node neighbor of our node
							bool good_neighbor = false;
							for (int t = 0; t<geom.node2face_count[n]; ++t)
								if (fneigh == geom.node2face_data[geom.node2face_start[n]+t])
								{
									good_neighbor = true;
									break;
								}

							if (good_neighbor)
							{
								// Compute edge ID
								const auto& edgeID = spade::utils::pairing_function(pf, fneigh);

								// Is this edge marked as a feature edge?
								const auto ipt = spade::algs::binary_search(edge_list, edgeID, invalid);
								if (ipt != invalid)
								{ // This is a feature edge. Add another dual vertex at midpoint
								
									// Dual vertex ID
									const auto inode = ipt + geom.faces.size(); // Counted after face centroids so increment the ID
										
									voronoi.face2node_data[voronoi.face2node_start[npoly]+nnode] = inode;
									++nnode;
								}
							}
						}
					}

					// If this is a vertex_node, one last add of ourselves
					if (map.type == vertex_node)
					{
						// Boundary intersection node iD
						const auto ipt   = spade::algs::binary_search(bndy_int, n, invalid);
						if (ipt == invalid)
						{
							print("Failed during dual connectivity formation when adding bndy_int node to list!");
							std::cin.get();
						}
								
						// Dual vertex ID
						const auto inode = ipt + edge_list.size() + geom.faces.size(); // Counted after face centroids and feature edges so increment the ID
								
						voronoi.face2node_data[voronoi.face2node_start[npoly]+nnode] = inode;
						++nnode;
					}
					
					// Re-order voronoi node ordering. Needs to be sequential and match nvec orientation of original surface
					reorder_voronoi_nodes(npoly, voronoi, geom);

					// Correct voronoi normal
					correct_voronoi_normal(npoly, voronoi, map.list[idual][0], geom);
					
					// Incremenent poly faces
					++npoly;
					
				} // End loop over primal vertex dual face list

			} // Primal vertex type check
			
		} // End primal vertex loop

		// Restricted Lloyd's Algorithm (Centroidal Voronoi Tesselation)
		restricted_lloyds_algorithm(nsmooth, geom, node_map, voronoi);
		
		// Complete
		print("Done!");
		return voronoi;
	}
}
