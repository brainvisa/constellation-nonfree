#!/usr/bin/env python
###############################################################################
# This software and supporting documentation are distributed by CEA/NeuroSpin,
# Batiment 145, 91191 Gif-sur-Yvette cedex, France. This software is governed
# by the CeCILL license version 2 under French law and abiding by the rules of
# distribution of free software. You can  use, modify and/or redistribute the
# software under the terms of the CeCILL license version 2 as circulated by
# CEA, CNRS and INRIA at the following URL "http://www.cecill.info".
###############################################################################

"""
This script does the following:
* 4-coloring maps

Main dependencies:

Author: Sandrine Lefranc, 2016
"""

# ---------------------------Imports-------------------------------------------


# System import
from __future__ import absolute_import
import random
from six.moves import range


# ---------------------------Functions-----------------------------------------

class Map:
    """
    """
    States = None
    Value = None

    def __init__(self, states, value):
        self.States = states
        self.Value = value


class Rule:
    """Connects two states indicating that they are adjacent.
    """
    Item = None
    Other = None
    Stringified = None

    def __init__(self, item, other, stringified):
        self.Item = item
        self.Other = other
        self.Stringified = stringified

    def __eq__(self, another):
        return hasattr(another, 'Item') and hasattr(another, 'Other') and \
            self.Item == another.Item and self.Other == another.Other

    def __hash__(self):
        return hash(self.Item) ^ hash(self.Other)

    def __str__(self):
        return self.Stringified


def generate_child(parent, colors, rules):
    """

    Parameters
    ----------
    parent: instance (mandatory)
    colors: list (mandatory)
    rules: list (mandatory)

    Returns
    -------
    Map(): instance
    """
    child_states = parent.States[:]
    index = random.randint(0, len(parent.States) - 1)
    state_id = random.randint(0, len(colors) - 1)
    child_states[index] = colors[state_id]
    value = get_value(child_states, rules)
    return Map(child_states, value)


def generate_parent(length, colors, rules):
    """

    Parameters
    ----------
    length: int (mandatory)
    colors: list (mandatory)
    rules: list (mandatory)

    Returns
    -------
    Map(): instance
    """
    child_states = []
    for i in range(0, length):
        state_id = random.randint(0, len(colors) - 1)
        child_states.append(colors[state_id])
    value = get_value(child_states, rules)
    return Map(child_states, value)


def get_best(rules, length, opt_value, colors):
    """

    Since the expected optimal situation will be that all adjacent states have
    different colors we can set the optimal value to the number of rules.

    Parameters
    ----------
    rules: list (mandatory)
    length: int (mandatory)
    opt_value: int (mandatory)
    colors: list (mandatory)

    Returns
    -------
    best_parent: instance
    """
    random.seed()
    best_parent = generate_parent(length, colors, rules)
    while best_parent.Value < opt_value:
        parent = generate_parent(length, colors, rules)
        attempts = 0
        while attempts < length:
            child = generate_child(parent, colors, rules)
            if child.Value > parent.Value:
                parent = child
                attempts = 0
            attempts += 1

        if best_parent.Value < parent.Value:
            best_parent, parent = parent, best_parent

    return best_parent


def get_value(colorlist, rules):
    """Returns an integer value representing how close that particular
    candidate comes to the optimal solution.

    Higher values are better.

    Parameters
    ----------
    colorlist: list (mandatory)
    rules: list (mandatory)

    Returns
    -------
    nb_selected_rules: int
    """
    nb_selected_rules = 0
    for rule in rules:
        if colorlist[rule.Item] != colorlist[rule.Other]:
            nb_selected_rules += 1

    return nb_selected_rules


def build_rules(items):
    """Build the set of rules.

    Whenever a state says it is adjacent to another state, the adjacent state
    should also say it is adjacent to the first state.

    Parameters
    ----------
    items: dict (mandatory)

    Returns
    -------
    rules: list
    """
    keys = sorted(items.keys())

    item2index = {}
    for idx, key in enumerate(keys):
        item2index[key] = idx

    add_rules = {}
    for key in keys:
        key_id = item2index[key]
        # get its neighbors
        adj_keys = items[key]
        for adj_key in adj_keys:
            adj_id = item2index[adj_key]
            temp = key_id
            if adj_id < temp:
                temp, adj_id = adj_id, temp
            rule_key = str(keys[temp]) + "->" + str(keys[adj_id])
            rule = Rule(temp, adj_id, rule_key)
            if rule in add_rules:
                add_rules[rule] += 1
            else:
                add_rules[rule] = 1

    rules = list(add_rules.keys())

    for k, v in add_rules.items():
        if v == 1:
            raise ValueError("rule %s is not bidirectional" % k)

    return rules
