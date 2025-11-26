from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional

from synplan.chem.reaction_routes.visualisation import cgr_display
from synplan.utils.visualisation import (
    get_route_svg_from_json,
    routes_clustering_report,
    routes_subclustering_report,
)

def compose_all_route_cgrs(tree_or_routes, route_id=None):
    """
    Process routes (reassign atom mappings) to compose RouteCGR.

    Parameters
    ----------
    tree_or_routes : synplan.mcts.tree.Tree
        or dict mapping route_id -> {step_id: ReactionContainer}
    route_id : int or None
        if None, do *all* winning nodes (or all keys of the dict);
        otherwise only that specific route.

    Returns
    -------
    dict or None
      - if route_id is None: {route_id: CGR, …}
      - if route_id is given: {route_id: CGR}
      - returns None on error
    """
    # dict-based branch
    if isinstance(tree_or_routes, dict):
        routes_dict = tree_or_routes

        def _single(route_id):
            res = compose_route_cgr(routes_dict, route_id)
            return res["cgr"] if res else None

        if route_id is not None:
            if route_id not in routes_dict:
                raise KeyError(f"Route {route_id} not in provided dict.")
            return {route_id: _single(route_id)}

        # all routes
        result = {route_id: _single(route_id) for route_id in sorted(routes_dict)}
        return result

    # tree-based branch
    tree = tree_or_routes
    route_cgrs = {}

    if route_id is not None:
        res = compose_route_cgr(tree, route_id)
        if res:
            route_cgrs[route_id] = res["cgr"]
        else:
            return None
        return route_cgrs

    for route_id in sorted(set(tree.winning_nodes)):
        res = compose_route_cgr(tree, route_id)
        if res:
            route_cgrs[route_id] = res["cgr"]

    return route_cgrs

def subcluster_all_clusters(groups, sb_cgrs_dict, route_cgrs_dict):
    """
    Subdivide each reaction cluster into detailed synthon-based subgroups.

    Iterates over all clusters in `groups`, applies `subcluster_one_cluster`
    to generate per-cluster synthons, then organizes routes by synthon detail.

    Parameters
    ----------
    groups : dict
        Mapping of cluster indices to cluster data.
    sb_cgrs_dict : dict
        Dictionary of SB-CGRs
    route_cgrs_dict : dict
        Dictionary of RoteCGRs

    Returns
    -------
    dict or None
        A dict mapping each cluster index to its subgroups dict,
        or None if any cluster fails to subcluster.
    """
    all_subgroups = {}
    for group_index, group in groups.items():
        group_synthons = subcluster_one_cluster(group, sb_cgrs_dict, route_cgrs_dict)
        if group_synthons is None:
            return None
        all_subgroups[group_index] = group_routes_by_synthon_detail(group_synthons)
    return all_subgroups

def cluster_routes(sb_cgrs: dict, use_strat=False):
    """
    Cluster routes objects based on their strategic bonds
      or CGRContainer object signature (not avoid mapping)

    Args:
        sb_cgrs: Dictionary mapping route_id to sb_cgr objects.

    Returns:
        Dictionary with groups keyed by '{length}.{index}' containing
        'sb_cgr', 'route_ids', and 'strat_bonds'.
    """
    temp_groups = defaultdict(
        lambda: {"route_ids": [], "sb_cgr": None, "strat_bonds": None}
    )

    # 1. Initial grouping based on the content of strategic bonds
    for route_id, sb_cgr in sb_cgrs.items():
        strat_bonds_list = extract_strat_bonds(sb_cgr)
        if use_strat == True:
            group_key = tuple(strat_bonds_list)
        else:
            group_key = str(sb_cgr)

        if not temp_groups[group_key]["route_ids"]:  # First time seeing this group
            temp_groups[group_key][
                "sb_cgr"
            ] = sb_cgr  # Store the first CGR as representative
            temp_groups[group_key][
                "strat_bonds"
            ] = strat_bonds_list  # Store the actual list

        temp_groups[group_key]["route_ids"].append(route_id)
        temp_groups[group_key][
            "route_ids"
        ].sort()  # Keep route_ids sorted for consistency

    for group_key in temp_groups.keys():
        temp_groups[group_key]["group_size"] = len(temp_groups[group_key]["route_ids"])

    # 2. Format the output dictionary with desired keys '{length}.{index}'
    final_grouped_results = {}
    group_indices = defaultdict(int)  # To track index for each length

    # Sort items by length of bonds first, then potentially by bonds themselves for consistent indexing
    # Sorting by the group_key (tuple of tuples) provides a deterministic order
    sorted_groups = sorted(
        temp_groups.items(), key=lambda item: (len(item[0]), item[0])
    )

    for group_key, group_data in sorted_groups:
        num_bonds = len(group_data["strat_bonds"])
        group_indices[num_bonds] += 1  # Increment index for this length (1-based)
        final_key = f"{num_bonds}.{group_indices[num_bonds]}"
        final_grouped_results[final_key] = group_data

    return final_grouped_results

def post_process_subgroup(
    subgroup,
):  # Under development: Error in replace_leaving_groups_in_synthon , 'cuz synthon_reaction.clean2d crashes
    """
    Drop leaving-groups common to all pathways and rebuild a minimal synthon.

    Scans the subgroup for leaving-groups present in every route, removes those
    from the CGR, re-assembles a clean ReactionContainer with the original core,
    updates `routes_data`, and flags the dict as processed.

    Parameters
    ----------
    subgroup : dict
        Must include keys for `routes_data` and the helpers
        (`all_lg_collect`, `find_const_lg`, etc.). If already
        post_processed, returns immediately.

    Returns
    -------
    dict
        The same dict, now with:
        - `'synthon_reaction'`: cleaned ReactionContainer
        - `'routes_data'`: filtered route table
        - `'post_processed'`: True
    """
    if "post_processed" in subgroup.keys() and subgroup["post_processed"] == True:
        return subgroup
    result = all_lg_collect(subgroup)
    # to find constant lg that need to be removed
    to_remove = [ind for ind, cgr_set in result.items() if len(cgr_set) == 1]
    new_synthon_cgr, new_lgs = replace_leaving_groups_in_synthon(subgroup, to_remove)
    synthon_reaction = ReactionContainer.from_cgr(new_synthon_cgr)
    synthon_reaction.clean2d()
    old_reactants = ReactionContainer.from_cgr(new_synthon_cgr).reactants
    target_mol = synthon_reaction.products[0]  # TO DO: target_mol might be non 0
    max_in_target_mol = max(target_mol._atoms)
    new_reactants = new_lg_reaction_replacer(
        synthon_reaction, new_lgs, max_in_target_mol
    )
    new_synthon_reaction = ReactionContainer(
        reactants=new_reactants, products=[target_mol]
    )
    new_synthon_reaction.clean2d()
    subgroup["synthon_reaction"] = new_synthon_reaction
    subgroup["routes_data"] = remove_and_shift(subgroup["routes_data"], to_remove)
    subgroup["post_processed"] = True
    subgroup["group_lgs"] = group_by_identical_values(subgroup["routes_data"])
    return subgroup


# ========== Jupyter / Colab 시각화 ==========
from IPython.display import display, HTML, SVG


# ========== Reduced CGR 유틸 함수들 (노트북에서 쓰던 것 그대로) ==========

def reduce_route_cgr_to_main_component(route_cgr):
    """
    하나의 Route CGR에서 '가장 큰 connected component'만 남긴 reduced CGR을 생성한다.

    Parameters
    ----------
    route_cgr
        CGR 객체 (예: CGRTools CGRContainer)
        compose_all_route_cgrs 의 value 로 나오는 객체.

    Returns
    -------
    reduced_cgr
        CGR 객체.
        - connected_components 가 여러 개인 경우: 가장 큰 component만 남긴 sub-CGR
        - connected_components 가 없거나 비어 있는 경우: 원본 route_cgr 그대로 반환
    """
    # SynPlanner + CGRTools 기준으로, connected_components 속성이 없을 가능성은 거의 없지만
    # 방어적으로 getattr 사용
    components = getattr(route_cgr, "connected_components", None)

    # component 정보가 없거나 비어 있으면 그냥 원본 반환
    if not components:
        return route_cgr

    # 각 component 는 "원자 인덱스 집합"이라고 가정하고,
    # 가장 원자 수가 많은 component를 route 대표 skeleton 으로 선택
    main_component = max(components, key=len)

    # 해당 component만 남긴 sub-CGR 생성
    reduced_cgr = route_cgr.substructure(main_component)

    return reduced_cgr


def compose_all_reduced_route_cgrs(source: Any, route_id: int | str | None = None) -> Dict[Any, Any]:
    """
    여러 Route CGR에 대해 'Reduced route CGR'을 일괄 생성한다.

    네 상황을 모두 커버하도록 설계:
    1) 이미 {route_id: CGR} 딕셔너리를 가진 경우
       - 예: all_route_cgrs = compose_all_route_cgrs(routes_dict_2)
       - compose_all_reduced_route_cgrs(all_route_cgrs)

    2) Tree 객체에서 바로 받고 싶은 경우
       - 예: compose_all_reduced_route_cgrs(tree)

    Parameters
    ----------
    source
        dict 또는 synplan.mcts.tree.Tree
        - dict: {route_id: CGR} 형태 (compose_all_route_cgrs 결과)
        - Tree: planning이 끝난 synplan.mcts.tree.Tree 객체

    route_id
        특정 route 하나만 대상으로 reduced CGR을 얻고 싶을 때 사용.
        None 이면 모든 route 에 대해 reduced CGR 생성.

    Returns
    -------
    Dict[route_id, reduced_CGR]
        {route_id: reduced_CGR} 형태의 딕셔너리.
        - route_id 를 지정한 경우: 해당 route_id 키만 있는 0 또는 1개짜리 dict
        - route_id 가 None 인 경우: 모든 route에 대한 reduced CGR dict
    """
    # 1) 이미 {route_id: CGR} dict 를 받은 경우
    if isinstance(source, dict):
        route_cgrs = source

    else:
        # 2) Tree 를 받은 경우: 먼저 route 별 CGR 생성
        #    (공식 SynPlanner route_cgr 모듈의 compose_all_route_cgrs 사용)
        from synplan.chem.reaction_routes.route_cgr import compose_all_route_cgrs

        route_cgrs = compose_all_route_cgrs(source)

    # 내부 헬퍼: 단일 route 하나를 reduced 로 변환
    def _reduce_single(rid):
        cgr = route_cgrs.get(rid)
        if cgr is None:
            return None
        return reduce_route_cgr_to_main_component(cgr)

    # 특정 route_id 만 요청한 경우
    if route_id is not None:
        if route_id not in route_cgrs:
            raise KeyError(f"Route {route_id} not found in provided CGR source.")
        reduced = _reduce_single(route_id)
        return {route_id: reduced} if reduced is not None else {}

    # 전체 route 에 대해 reduced CGR 생성
    reduced_dict: Dict[Any, Any] = {}

    # 정렬 가능한 경우에는 정렬해서, 아니라면 그냥 순회
    try:
        route_ids = sorted(route_cgrs)
    except TypeError:
        route_ids = list(route_cgrs.keys())

    for rid in route_ids:
        reduced = _reduce_single(rid)
        if reduced is not None:
            reduced_dict[rid] = reduced

    return reduced_dict

