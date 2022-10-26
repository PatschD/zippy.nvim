# zippy.nvim

nvim plugin written in lua to help with debugging. this is inspired by the popular vscode extension [turbo console log](https://marketplace.visualstudio.com/items?itemName=ChakrounAnas.turbo-console-log), for javascript, lua and python.

the main functionality is to insert meaningful print statements automatically.

![image](https://user-images.githubusercontent.com/80820813/198012225-9cafe569-fa80-461e-9a58-dbbb25653f50.png)

![image](https://user-images.githubusercontent.com/80820813/198012504-11da9e43-4bdd-4982-b37f-7c8670e91d37.png)

the print statement is formatted as 
```
- file name -> line number -> variable
```

if you have [nvim navic](https://github.com/SmiteshP/nvim-navic) installed, breadcrumbs to your function are additionally added:
```
- file name -> line number -> breadcrumbs -> variable
```

## Installation

Packer:

```
use 'PatschD/zippy.nvim'
```

## Setup

<a href="https://github.com/nvim-treesitter/nvim-treesitter" target="_blank" rel="noopener noreferrer"> Tree-sitter is required</a> <br/>

<b> if you are happy with the defaults, nothing has to be configured.</b>

otherwise two options can be configured:

- formatting:

as lua, javascript and python are quite different (indentations, newlines), zippy.nvim focuses on creating valid code, then calls the lsp to format. this can be disabled, performed only on the affected rows or perfomed globally. Not all formatters can format only a certain range (such as black for python)

- print function:
if you have a custom print function you can pass that instead.

```
require("zippy").setup({
  ["lua"] = {
    language_format_behaviour = "none", 
    print_function = "print",
  },
  ["typescriptreact] = {
    print_function = "console.warn",
  }
)
```

## Keymap

vim.keymap.set("n", "<leader>lg", "<cmd>lua require('zippy').insert_print()<CR>")


