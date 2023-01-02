import flet as ft
from utilities import *
from flet import theme

def main(page: ft.Page):
    page.title = "Sequence Alignment"
    page.vertical_alignment = ft.MainAxisAlignment.SPACE_EVENLY
   

    #page.vertical_alignment = ft.MainAxisAlignment.SPACE_AROUND

    Sequence1 = ft.TextField(label="Enter Sequence 1", autofocus=True)
    Sequence2 = ft.TextField(label="Enter Sequence 2", autofocus=True)

    Sequences = ft.Column()
    cg = ft.RadioGroup(content=ft.Row([
    ft.Radio(value="DNA", label="DNA"),
    ft.Radio(value="RNA", label="RNA"),
    ft.Radio(value="PROTEIN", label="PROTEIN")]))
    match = ft.TextField(label="MATCH",value="0", text_align=ft.TextAlign.CENTER, width=100,autofocus=True)
    mismatch = ft.TextField(label="MISMATCH",value="0", text_align=ft.TextAlign.CENTER, width=120,autofocus=True)
    gap = ft.TextField(label="GAP",value="0", text_align=ft.TextAlign.CENTER, width=100)


    selected_files = ft.Text()

    def pick_files_result(e: ft.FilePickerResultEvent):
        selected_files.value = (
            ", ".join(map(lambda f: f.name, e.files)) if e.files else "Cancelled!"
        )
        selected_files.update()

    pick_files_dialog = ft.FilePicker(on_result=pick_files_result)

    page.overlay.append(pick_files_dialog)

    def btn_click(e):
        Sequences.controls.clear()
        result1, result2 = global_alignment(Sequence1.value, Sequence2.value, match.value, mismatch.value, gap.value)

        Sequences.controls.append(ft.Text(f"Sequence 1 :  {result1} "))
        Sequences.controls.append(ft.Text(f"Sequence 2 :  {result2} "))
        page.update()


   # MATCH COUNT
    def match_minus_click(e):
        match.value = str(int(match.value) - 1)
        page.update()

    def match_plus_click(e):
        match.value = str(int(match.value) + 1)
        page.update()
    

    # MISMATCH COUNT
    def mismatch_minus_click(e):
        mismatch.value = str(int(mismatch.value) - 1)
        page.update()

    def mismatch_plus_click(e):
        mismatch.value = str(int(mismatch.value) + 1)
        page.update()
    
    # GAP COUNT
    def gap_minus_click(e):
        gap.value = str(int(gap.value) - 1)
        page.update()

    def gap_plus_click(e):
        gap.value = str(int(gap.value) + 1)
        page.update()
        
    page.theme = theme.Theme(color_scheme_seed="red")
    page.update()


    page.add(
        ft.Row(
            [
                ft.ElevatedButton(
                    "Pick files",
                    icon=ft.icons.UPLOAD_FILE,
                    on_click=lambda _: pick_files_dialog.pick_files(
                        allow_multiple=True)),
                          selected_files,
            ]
        ),
        cg,
        ft.Row(
            [
                ft.IconButton(ft.icons.REMOVE, on_click=match_minus_click),
                match,
                ft.IconButton(ft.icons.ADD, on_click=match_plus_click),


                 ft.IconButton(ft.icons.REMOVE, on_click=mismatch_minus_click),
                mismatch,
                ft.IconButton(ft.icons.ADD, on_click=mismatch_plus_click),



                  ft.IconButton(ft.icons.REMOVE, on_click=gap_minus_click),
                gap,
                ft.IconButton(ft.icons.ADD, on_click=gap_plus_click)
            ]
            
        ),
        ft.Column(
        [
        Sequence1,
        Sequence2,
        ft.Row(

        [
           ft.ElevatedButton("Global Alignment", on_click=btn_click), 
           ft.ElevatedButton("Local Alignment", on_click=btn_click),
        ]
        ) 
        ,
        Sequences])
    )
    

ft.app(target=main)
#ft.app(target=main, view=ft.WEB_BROWSER)

