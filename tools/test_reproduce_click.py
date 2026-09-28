import sys
from acu3d.model import Atlas
from acu3d.surface_viewer import SurfaceViewer

def test():
    print("Initializing Atlas & SurfaceViewer...")
    atlas = Atlas()
    sv = SurfaceViewer(atlas, off_screen=True)
    
    print("\n--- Test cycle_yuan_luo ---")
    sv.cycle_yuan_luo()
    print("Selected point:", sv.selected_point)
    print("Highlight actor pair vis:", sv.highlight_actor_pair.GetVisibility() if sv.highlight_actor_pair else None)
    print("Relation tube actor vis:", sv.relation_tube_actor.GetVisibility() if sv.relation_tube_actor else None)





if __name__ == "__main__":
    test()
