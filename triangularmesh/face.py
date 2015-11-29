class Face(object):
    def __init__(self, node_list, diameter = 0, parent_diameter = 0, parent_index = -1):
        if type(node_list) is not int:
            self.nodes = node_list
            self.parent_diameter = parent_diameter
            self.diameter = diameter
        else:
            self.nodes = [node_list, diameter, parent_diameter]
            self.parent_diameter = 0
            self.diameter = 0

        self.index = -1
        self.parent_index = parent_index

    def __eq__(self, other):
        return set(self) == set(other)

    def __str__(self):
        return '[' + str(self.nodes[0]) + ',' + str(self.nodes[1]) + ',' + str(self.nodes[2]) + ']'

    def __repr__(self):
        return '[' + str(self.nodes[0]) + ',' + str(self.nodes[1]) + ',' + str(self.nodes[2]) + ']'

    def __getitem__(self, item):
        return self.nodes[item]

    def __setitem__(self, key, value):
        self.nodes[key] = value

    class __metaclass__(type):
        def __iter__(self):
            for t in self.nodes:
                yield t
