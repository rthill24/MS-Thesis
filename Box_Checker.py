# this module checks if the stiffened box stiffeners intersect and returns ABS factors
# author: Richard Thill
# date: 10-29-2024

import TPanel_trans

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

        if self.nstiff_side == 1:
            stiff_spacing = self.B_side/2
        else:
            stiff_spacing = (self.B_side)/(self.nstiff_side-1)
        
        h_bot_side_stiff = self.tp_bot + stiff_spacing
        h_bot_side_stiff_flange = h_bot_side_stiff - (self.tw_side/2) - (self.bf_side/2)
        h_bot_stiff = self.tp_bot + self.hw_top + self.tf_top
        #print (stiff_spacing, h_bot_side_stiff, h_bot_side_stiff_flange, h_bot_stiff) 
        if h_bot_stiff >= h_bot_side_stiff_flange:
            return bool(1) 
        else:
            return bool(0)

    def check_geometry(self):
        valid_bot_data = TPanel_trans.TPanel_trans(self.test_panel_bot)
        valid_bot = valid_bot_data.geoValid()

        valid_side_data = TPanel_trans.TPanel_trans(self.test_panel_side)
        valid_side = valid_side_data.geoValid()

        valid_top_data = TPanel_trans.TPanel_trans(self.test_panel_top)
        valid_top = valid_top_data.geoValid()
        #print(valid_bot, valid_side, valid_top)