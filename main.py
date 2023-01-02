import flet as ft
from utilities import *

def main(page: ft.Page):
    page.title = "Sequence Alignment"
    #page.vertical_alignment = ft.MainAxisAlignment.SPACE_AROUND

    Sequence1 = ft.TextField(label="Enter Sequence1", autofocus=True)
    Sequence2 = ft.TextField(label="Enter Sequence2", autofocus=True)

    Sequences = ft.Column()
    DNA = ft.Radio(label="DNA", value=False)
    RNA = ft.Radio(label="RNA", value=False)
    PROTEIN = ft.Radio(label="PROTEIN", value=False)
    match_text=ft.Text("Match")
    dismatch_text=ft.Text("Dismatch")
    gap_text=ft.Text("Gap")

    match = ft.TextField(value="0", text_align=ft.TextAlign.RIGHT, width=50)
    dismatch = ft.TextField(value="0", text_align=ft.TextAlign.RIGHT, width=50)
    gap = ft.TextField(value="0", text_align=ft.TextAlign.RIGHT, width=50)


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
        Sequences.controls.append(ft.Text(f"Sequence 1 are :  {Sequence1.value} "))
        Sequences.controls.append(ft.Text(f"Sequence 2 are :  {Sequence2.value} "))

        result1, result2 = global_alignment(Sequence1.value, Sequence2.value)

        Sequences.controls.append(ft.Text(f"Sequence 1 are :  {result1} "))
        Sequences.controls.append(ft.Text(f"Sequence 2 are :  {result2} "))
        page.update()


   # MATCH COUNT
    def match_minus_click(e):
        match.value = str(int(match.value) - 1)
        page.update()

    def match_plus_click(e):
        match.value = str(int(match.value) + 1)
        page.update()
    

    # DISMATCH COUNT

    def dismatch_minus_click(e):
        dismatch.value = str(int(dismatch.value) - 1)
        page.update()

    def dismatch_plus_click(e):
        dismatch.value = str(int(dismatch.value) + 1)
        page.update()
    
    # GAP COUNT
    def gap_minus_click(e):
        gap.value = str(int(gap.value) - 1)
        page.update()

    def gap_plus_click(e):
        gap.value = str(int(gap.value) + 1)
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
        ft.Row([DNA,RNA,PROTEIN]),
        ft.Row(
            [
                match_text,ft.IconButton(ft.icons.REMOVE, on_click=match_minus_click),
                match,
                ft.IconButton(ft.icons.ADD, on_click=match_plus_click),


                 dismatch_text,ft.IconButton(ft.icons.REMOVE, on_click=dismatch_minus_click),
                dismatch,
                ft.IconButton(ft.icons.ADD, on_click=dismatch_plus_click),



                  gap_text,ft.IconButton(ft.icons.REMOVE, on_click=gap_minus_click),
                gap,
                ft.IconButton(ft.icons.ADD, on_click=gap_plus_click)

                
            ]  ,alignment=ft.MainAxisAlignment.CENTER,
            
        ),
        ft.Column(
        [
        Sequence1,
        Sequence2,
        ft.ElevatedButton("show Alignment", on_click=btn_click),
        Sequences]
        )
    )
  
ft.app(target=main)
#ft.app(target=main, view=ft.WEB_BROWSER)

