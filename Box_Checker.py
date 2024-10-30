# this module checks if the stiffened box stiffeners intersect and returns ABS factors
# author: Richard Thill
# date: 10-29-2024

class boxchecker:

    def __init__(self, tp_bot, hw_top, tf_top, B_side, nstiff_side, tw_side, bf_side, test_panel_bot, test_panel_side, test_panel_top):
        self.tp_bot = tp_bot
        self.hw_top = hw_top
        self.tf_top = tf_top
        self.B_side = B_side
        self.nstiff_side = nstiff_side
        self.tw_side = tw_side
        self.bf_side = bf_side
        self.test_panel_bot = test_panel_bot
        self.test_panel_side = test_panel_side
        self.test_panel_top = test_panel_top

    def check_corner_stiffs(self):
        stiff_spacing = (self.B_side)/(self.nstiff_side-1)
        h_bot_side_stiff = self.tp_bot + stiff_spacing
        h_bot_side_stiff_flange = h_bot_side_stiff - (self.tw_side/2) - (self.bf_side/2)
        h_bot_stiff = self.tp_bot + self.hw_top + self.tf_top
        print (stiff_spacing, h_bot_side_stiff, h_bot_side_stiff_flange, h_bot_stiff) 
        if h_bot_stiff >= h_bot_side_stiff_flange:
            print ("Corner stiffeners are intersecting") 
        else:
            print ("Corner stiffeners are not intersecting")

    def check_geometry(self):
        valid_bot = self.test_panel_bot.geoValid()
        valid_side = self.test_panel_side.geoValid()
        valid_top = self.test_panel_top.geoValid()
        print(valid_bot, valid_side, valid_top)